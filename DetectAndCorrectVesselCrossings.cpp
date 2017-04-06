#include "mex.h"

#include <exception>
#include <set>
#include <numeric>
#include <cassert>

#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkAdjacentVertexIterator.h>
#include <vtkEdgeListIterator.h>
#include <vtkOutEdgeIterator.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkGraphToPolyData.h>

class SearchRegion
{
public:
    SearchRegion(double* originCoords, double radius):
        originCoords(originCoords), radiusSquared(radius*radius)
    {
    }

    bool IsPointContained(double* pointCoords) const
    {
        double distanceSquared = vtkMath::Distance2BetweenPoints(pointCoords, originCoords);            
        return(distanceSquared <= radiusSquared);
    }

    double GetDistanceSquaredInterpolator(double* pointCoords) const
    {
        assert(IsPointContained(pointCoords));        
        double distanceSquared = vtkMath::Distance2BetweenPoints(pointCoords, originCoords);            
        return (1 - distanceSquared / radiusSquared);
    }

private:
    double* originCoords;
    double radiusSquared;
};


class NetworkAnalyser
{
public:
    NetworkAnalyser(const std::string& vtpSourceFilename)
    {
        // Read VTK polydata from file (.vtp file)
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(vtpSourceFilename.c_str());
        reader->Update();
        networkSkeleton = reader->GetOutput();

        // Create VTK graph datastructure
        networkGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

        numVertices = networkSkeleton->GetNumberOfPoints();    
        for (unsigned vertex_id=0; vertex_id<numVertices; ++vertex_id)
        {
            networkGraph->AddVertex();
        }

        vtkCellArray* edges = networkSkeleton->GetLines();
        vtkSmartPointer<vtkIdList> edge = vtkSmartPointer<vtkIdList>::New();
        while(edges->GetNextCell(edge))
        {
            networkGraph->AddGraphEdge(edge->GetId(0), edge->GetId(1));
        }

        vertexCoordinates.resize(numVertices);
        vertexRadius.resize(numVertices);
        for (unsigned vertex_id=0; vertex_id<numVertices; ++vertex_id)
        {
            vertexCoordinates[vertex_id] = new double[3];
            networkSkeleton->GetPoints()->GetPoint(vertex_id, vertexCoordinates[vertex_id]);
            vertexRadius[vertex_id] = networkSkeleton->GetPointData()->GetScalars()->GetTuple(vertex_id)[0];
        }    
        
    }
    
    ~NetworkAnalyser()
    {
        for (std::vector<double*>::iterator iter=vertexCoordinates.begin();
             iter!=vertexCoordinates.end();
             ++iter)
        {
            delete[] *iter;
        }

    }
    
    void FindAndCorrectDegree4VesselCrossings()
    {
        // Loop over all the vertices and find those with degree 4
        for (unsigned vertex_id=0; vertex_id<numVertices; ++vertex_id)
        {
            if(networkGraph->GetDegree(vertex_id) == 4)
            {
                std::cout << "Degree 4 vessel crossing at vertex: " << vertex_id << std::endl;

                // Get radius and location of original branching point
                double radius = vertexRadius[vertex_id];
                double original_coords[3];
                original_coords[0] = vertexCoordinates[vertex_id][0];
                original_coords[1] = vertexCoordinates[vertex_id][1];
                original_coords[2] = vertexCoordinates[vertex_id][2];

                // Add vertex to the graph
                unsigned new_vertex = networkGraph->AddVertex();
                numVertices++;

                // Store its location
                double* new_vertex_coords = new double[3]; 
                new_vertex_coords[0] = vertexCoordinates[vertex_id][0];
                new_vertex_coords[1] = vertexCoordinates[vertex_id][1];
                new_vertex_coords[2] = vertexCoordinates[vertex_id][2];
                vertexCoordinates.push_back(new_vertex_coords);

                // Store its radius
                vertexRadius.push_back(radius);

                // Modify graph connectivity (add/remove edges to split the crossing into two segments that run at differetn depths)
                std::vector<unsigned> vertex_neighbours;
                std::vector<unsigned> connecting_edges;
                for (unsigned neigh_id=0; neigh_id<4; ++neigh_id)
                {
                    vtkOutEdgeType graph_edge = networkGraph->GetOutEdge(vertex_id, neigh_id);
                    vertex_neighbours.push_back(graph_edge.Target);
                    connecting_edges.push_back(graph_edge.Id);
                }

                try
                {
                    std::pair<unsigned, unsigned> aligned_neighbours_indices = GetPairOfAlignedVertices(vertex_id, vertex_neighbours);

                    networkGraph->RemoveEdge(connecting_edges[aligned_neighbours_indices.first]);
                    networkGraph->RemoveEdge(connecting_edges[aligned_neighbours_indices.second]);

                    networkGraph->AddEdge(vertex_neighbours[aligned_neighbours_indices.first], new_vertex);
                    networkGraph->AddEdge(vertex_neighbours[aligned_neighbours_indices.second], new_vertex);
                }
                catch (std::exception &e)
                {
                    std::cout << "Not possible to classify neighbours of vertex " << vertex_id << std::endl;
                    exit(-1);
                }

                // Compute mean radii of the new segment containing the new vertex 
                std::set<unsigned> visited_neighbours;
                std::vector<double> radii;
                GetRadiiWithinDistance(new_vertex, SearchRegion(original_coords, 2*radius), visited_neighbours, radii);
                double mean_radii= std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();

                // Adjust vertex location and radii for nodes along segment containing the new vertex 
                visited_neighbours.clear();
                double displacement = radius; // displacing by the original vertex radius seems to be enough to avoid kissing segments
                AdjustLocationAndRadiiOfNeighbouringVertices(new_vertex, displacement, mean_radii, 
                                                             SearchRegion(original_coords, 2*radius), visited_neighbours);

                // Compute mean radii of the new segment containing the original vertex 
                visited_neighbours.clear();
                radii.clear();
                GetRadiiWithinDistance(vertex_id, SearchRegion(original_coords, 2*radius), visited_neighbours, radii);
                mean_radii= std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();

                // Adjust vertex location and radii for nodes along segment containing the new vertex 
                visited_neighbours.clear();
                AdjustLocationAndRadiiOfNeighbouringVertices(vertex_id, -displacement, mean_radii, 
                                                             SearchRegion(original_coords, 2*radius), visited_neighbours);
            }
        }
    }

    typedef std::pair<unsigned, unsigned> UnsignedPair;
    std::pair<UnsignedPair, UnsignedPair> ComputeBranchContinuity(unsigned firstCrossing, unsigned secondCrossing)
    {
        vtkSmartPointer<vtkAdjacentVertexIterator> neighbour_iterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();

        std::set<unsigned> first_crossing_neighbours;
        networkGraph->GetAdjacentVertices(firstCrossing, neighbour_iterator);
        while (neighbour_iterator->HasNext())
        {
            first_crossing_neighbours.insert(neighbour_iterator->Next());
        }

        std::set<unsigned> second_crossing_neighbours;
        networkGraph->GetAdjacentVertices(secondCrossing, neighbour_iterator);
        while (neighbour_iterator->HasNext())
        {
            second_crossing_neighbours.insert(neighbour_iterator->Next());
        }
        
        std::vector<unsigned> path_first_to_second(ComputePathBetweenTwoVertices(firstCrossing, secondCrossing));
        for (auto node_in_path_iter=path_first_to_second.begin(); node_in_path_iter!=path_first_to_second.end(); ++node_in_path_iter)
        {
            first_crossing_neighbours.erase(*node_in_path_iter);
            second_crossing_neighbours.erase(*node_in_path_iter);
        }

        assert(first_crossing_neighbours.size() == 2);
        assert(second_crossing_neighbours.size() == 2);        
        
        auto first_crossing_neighbours_iter = first_crossing_neighbours.begin();
        unsigned point_a = *first_crossing_neighbours_iter;
        ++first_crossing_neighbours_iter;
        unsigned point_b = *first_crossing_neighbours_iter;        

        auto second_crossing_neighbours_iter = second_crossing_neighbours.begin();
        unsigned point_c = *second_crossing_neighbours_iter;
        ++second_crossing_neighbours_iter;
        unsigned point_d = *second_crossing_neighbours_iter;                
        
        std::pair<UnsignedPair, UnsignedPair> paired_up_nodes;
        
        paired_up_nodes.first.first = point_a;
        paired_up_nodes.second.first = point_b;

        double* point_first_coords = vertexCoordinates[firstCrossing];
        double* point_a_coords = vertexCoordinates[point_a];
        double* point_b_coords = vertexCoordinates[point_b];

        double* point_second_coords = vertexCoordinates[secondCrossing];
        double* point_c_coords = vertexCoordinates[point_c];
        double* point_d_coords = vertexCoordinates[point_d];
        
        double vector_a_edge[3];
        vtkMath::Subtract(point_a_coords, point_first_coords, vector_a_edge);
        vtkMath::Normalize(vector_a_edge);

        double vector_c_edge[3];
        vtkMath::Subtract(point_c_coords, point_second_coords, vector_c_edge);
        vtkMath::Normalize(vector_c_edge);
        
        double vector_d_edge[3];
        vtkMath::Subtract(point_d_coords, point_second_coords, vector_d_edge);
        vtkMath::Normalize(vector_d_edge);
        
        double a_edge_dot_c_edge = fabs(vtkMath::Dot(vector_a_edge, vector_c_edge));
        double a_edge_dot_d_edge = fabs(vtkMath::Dot(vector_a_edge, vector_d_edge));
        
        if (a_edge_dot_c_edge > a_edge_dot_d_edge)
        {
            paired_up_nodes.first.second = point_c;
            paired_up_nodes.second.second = point_d;
        }
        else
        {
            paired_up_nodes.first.second = point_d;
            paired_up_nodes.second.second = point_c;
        }
        
        return paired_up_nodes;
    }
    
    void Correct2PairsOfDegree3VesselCrossings(std::vector<std::pair<unsigned, unsigned> > crossingList)
    {
        for (auto iter=crossingList.begin(); iter!=crossingList.end(); ++iter)
        {
            unsigned first_crossing = iter->first;
            unsigned second_crossing = iter->second;
            assert (first_crossing < numVertices);
            assert (second_crossing < numVertices);
            
            // Compute the pairs of neighbours of first_crossing and second_crossing that will be linked in each of the new branches
            std::pair<UnsignedPair, UnsignedPair> paired_up_nodes = ComputeBranchContinuity(first_crossing, second_crossing);
            
            // first_crossing_neigh_a will connect with second_crossing_neigh_a via the original path between first_crossing and seccon_crossing
            // first_crossing_neigh_b will connect with second_crossing_neigh_b via a new set of vertices to be created
            unsigned first_crossing_neigh_a = paired_up_nodes.first.first;
            unsigned first_crossing_neigh_b = paired_up_nodes.second.first;
            unsigned second_crossing_neigh_a = paired_up_nodes.first.second;
            unsigned second_crossing_neigh_b = paired_up_nodes.second.second;

            /*
             * Remove (first_crossing_neigh_b, first_crossing) edge
             */
            unsigned edge_id = GetEdgeForNodePair(first_crossing_neigh_b, first_crossing);
            networkGraph->RemoveEdge(edge_id);
            
            /*
             * Create the node that overlaps with first_crossing on the new branch
             * and connect it to first_crossing_neigh_b
             */
            unsigned new_first_crossing = networkGraph->AddVertex();
            numVertices++;

            // Store its location
            double* new_vertex_coords = new double[3]; 
            new_vertex_coords[0] = vertexCoordinates[first_crossing][0];
            new_vertex_coords[1] = vertexCoordinates[first_crossing][1];
            new_vertex_coords[2] = vertexCoordinates[first_crossing][2];
            vertexCoordinates.push_back(new_vertex_coords);

            // Store its radius
            vertexRadius.push_back(vertexRadius[first_crossing]);
            
            // Add edge
            networkGraph->AddEdge(first_crossing_neigh_b, new_first_crossing);
            
            /*
             * Create the new branch that will connect first_crossing_neigh_b and second_crossing_neigh_b
             */
            std::vector<unsigned> vertices_between_crossings = ComputePathBetweenTwoVertices(first_crossing, second_crossing);
            vertices_between_crossings[0] = new_first_crossing;
            
            for (unsigned vertex_id=1; vertex_id<vertices_between_crossings.size(); ++vertex_id)
            {
                unsigned source = vertices_between_crossings[vertex_id-1];
                unsigned target = vertices_between_crossings[vertex_id];
                
                // Create a node that overlaps target (the one overlaping source was already created)
                unsigned new_vertex = networkGraph->AddVertex();
                numVertices++;

                // Store its location
                double* new_vertex_coords = new double[3]; 
                new_vertex_coords[0] = vertexCoordinates[target][0];
                new_vertex_coords[1] = vertexCoordinates[target][1];
                new_vertex_coords[2] = vertexCoordinates[target][2];
                vertexCoordinates.push_back(new_vertex_coords);

                // Store its radius
                vertexRadius.push_back(vertexRadius[target]);

                // Add edge
                networkGraph->AddEdge(source, new_vertex);
                
                // Update data structure for next iteration
                vertices_between_crossings[vertex_id] = new_vertex;
            }

            /*
             * Remove (second_crossing, second_crossing_neigh_b) edge
             */
            edge_id = GetEdgeForNodePair(second_crossing, second_crossing_neigh_b);
            networkGraph->RemoveEdge(edge_id);
            
            /*
             * Connect second_crossing_neigh_b with the new branch
             */
            unsigned new_second_crossing = vertices_between_crossings[vertices_between_crossings.size()-1];
            networkGraph->AddEdge(new_second_crossing, second_crossing_neigh_b);

            // Compute mean radii of the new segment containing the new vertex 
            std::set<unsigned> visited_neighbours;
            std::vector<double> radii;
            double radius = vertexRadius[first_crossing];
            double original_coords[3];
            original_coords[0] = vertexCoordinates[first_crossing][0];
            original_coords[1] = vertexCoordinates[first_crossing][1];
            original_coords[2] = vertexCoordinates[first_crossing][2];
            GetRadiiWithinDistance(first_crossing, SearchRegion(original_coords, 2*radius), visited_neighbours, radii);
            double mean_radii= std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();

            // Adjust vertex location and radii for nodes along segment containing the new vertex 
            visited_neighbours.clear();
            double displacement = 1.5*radius; // displacing by the original vertex radius seems to be enough to avoid kissing segments
            AdjustLocationAndRadiiOfNeighbouringVertices(first_crossing, displacement, mean_radii, 
                                                         SearchRegion(original_coords, 5*radius), visited_neighbours);

            // Compute mean radii of the new segment containing the original vertex 
            visited_neighbours.clear();
            radii.clear();
            radius = vertexRadius[new_second_crossing];
            original_coords[0] = vertexCoordinates[new_second_crossing][0];
            original_coords[1] = vertexCoordinates[new_second_crossing][1];
            original_coords[2] = vertexCoordinates[new_second_crossing][2];            
            GetRadiiWithinDistance(new_second_crossing, SearchRegion(original_coords, 2*radius), visited_neighbours, radii);
            mean_radii= std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();

            // Adjust vertex location and radii for nodes along segment containing the new vertex 
            visited_neighbours.clear();
            displacement = 1.5*radius;
            AdjustLocationAndRadiiOfNeighbouringVertices(new_second_crossing, -displacement, mean_radii, 
                                                         SearchRegion(original_coords, 5*radius), visited_neighbours);            
        }        
    }
    
    unsigned GetEdgeForNodePair(unsigned nodeA, unsigned nodeB)
    {
        vtkSmartPointer<vtkOutEdgeIterator> edge_iter = vtkSmartPointer<vtkOutEdgeIterator>::New();
        networkGraph->GetOutEdges(nodeA, edge_iter);
        while(edge_iter->HasNext())
        {
            vtkOutEdgeType graph_edge = edge_iter->Next();
            if (graph_edge.Target == nodeB)
            {
                return graph_edge.Id;
            }
        }
        
        assert(false);
        return UINT_MAX;
    }
    
    std::vector<unsigned> ComputePathBetweenTwoVertices(unsigned sourceVertex, unsigned targetVertex)
    {
        // Use Dijkstra's algorithm to compute path between sourceVertex and targetVertex
        vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
        dijkstra->SetInputData(networkSkeleton);
        
        // The algorithm seems to returns the path from SetEndVertex to SetStartVertex
        dijkstra->SetStartVertex(targetVertex);
        dijkstra->SetEndVertex(sourceVertex);
        
        dijkstra->Update();
        vtkIdList* computed_path = dijkstra->GetIdList();

        std::vector<unsigned> path;
        for(unsigned node_index=0; node_index<computed_path->GetNumberOfIds(); ++node_index)
        {
            path.push_back(computed_path->GetId(node_index));
        }

        // Ensure path is from sourceVertex to targetVertex
        assert(path[0] == sourceVertex);
        assert(path[1] == targetVertex);
        
        return path;
    }    
    
    void WriteToDisc(const std::string& vtpOutputFilename)
    {
        // Create a polydata object
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();    

        // Create and fill a point container and add points to polydata
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (unsigned vertex=0; vertex<numVertices; vertex++)
        {
            points->InsertNextPoint(vertexCoordinates[vertex][0], 
                                    vertexCoordinates[vertex][1], 
                                    vertexCoordinates[vertex][2]);
        }
        polydata->SetPoints(points);

        // Allocate memory and add the definition of edges to the polydata object
        polydata->Allocate();

        vtkSmartPointer<vtkEdgeListIterator> edge_iterator = vtkSmartPointer<vtkEdgeListIterator>::New();
        networkGraph->GetEdges(edge_iterator);
        while(edge_iterator->HasNext())
        {
            vtkEdgeType edge = edge_iterator->Next();

            vtkIdType edge_ends[2];
            edge_ends[0] = edge.Source; 
            edge_ends[1] = edge.Target;

            polydata->InsertNextCell(VTK_LINE, 2, edge_ends);
        }

        // Create and fill the array containing the plexus radii at each vertex.
        vtkSmartPointer<vtkDoubleArray> radii_array = vtkSmartPointer<vtkDoubleArray>::New();
        radii_array->SetNumberOfComponents(1);
        radii_array->SetName("Radius");    
        for (unsigned vertex_id=0; vertex_id<numVertices; vertex_id++)
        {
            radii_array->InsertNextValue(vertexRadius[vertex_id]);
        }
        polydata->GetPointData()->AddArray(radii_array);
        polydata->GetPointData()->SetActiveScalars("Radius");

        WritePolyDataToDisk(polydata, vtpOutputFilename);        
    }
    
private:

    void WritePolyDataToDisk(const vtkSmartPointer<vtkPolyData> polyData, const std::string& fileName)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(polyData);

        // Optional - set the mode. The default is binary.
        //writer->SetDataModeToBinary();
        //writer->SetDataModeToAscii();

        writer->Write();
    }
    
    std::pair<unsigned, unsigned> GetPairOfAlignedVertices(unsigned vertex_id, std::vector<unsigned> vertex_neighbours)
    {
        // Take the first neighbour in the list and try to find another neighbour that is aligned with it.
        double point_a_coords[3];
        networkSkeleton->GetPoints()->GetPoint(vertex_neighbours[0], point_a_coords);

        double point_b_coords[3];
        networkSkeleton->GetPoints()->GetPoint(vertex_id, point_b_coords);

        double vector_ab[3];
        vtkMath::Subtract(point_b_coords, point_a_coords, vector_ab);
        vtkMath::Normalize(vector_ab);

        unsigned aligned_neighbour = UINT_MAX;
        for(unsigned neigh_id=1; neigh_id<vertex_neighbours.size(); ++neigh_id)
        {
            double point_c_coords[3];
            networkSkeleton->GetPoints()->GetPoint(vertex_neighbours[neigh_id], point_c_coords);

            double vector_bc[3];
            vtkMath::Subtract(point_c_coords, point_b_coords, vector_bc);
            vtkMath::Normalize(vector_bc);

            double ab_dot_bc = vtkMath::Dot(vector_ab, vector_bc);
            if (fabs(1.0-ab_dot_bc) < 1e-3)
            {
                aligned_neighbour = neigh_id;
                break;
            }        
        }
        if (aligned_neighbour == UINT_MAX)
        {
            throw std::exception();
        }

        return std::pair<unsigned, unsigned>(0, aligned_neighbour);
    }    

    void AdjustLocationAndRadiiOfNeighbouringVertices(unsigned vertexId, double zDisplacement, double newRadius, const SearchRegion& searchRegion, std::set<unsigned>& visitedNeighbours)
    {
        // Check if we have visited this vertex in the recursion
        if (visitedNeighbours.find(vertexId) != visitedNeighbours.end())
        {
            return;
        }
        // Mark vertex as visited, if not
        visitedNeighbours.insert(vertexId);
        
        // Displace in the z coordinate
        vertexCoordinates[vertexId][2] += zDisplacement * searchRegion.GetDistanceSquaredInterpolator(vertexCoordinates[vertexId]);
        vertexRadius[vertexId] = newRadius;
        
        // Recursively visit neighbouring nodes contained in the search region
        vtkSmartPointer<vtkAdjacentVertexIterator> neighbour_iterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        networkGraph->GetAdjacentVertices(vertexId, neighbour_iterator);
        while (neighbour_iterator->HasNext())
        {
            unsigned neighIndex = neighbour_iterator->Next();
            
            if (searchRegion.IsPointContained(vertexCoordinates[neighIndex]))
            {
                AdjustLocationAndRadiiOfNeighbouringVertices(neighIndex, zDisplacement, newRadius, searchRegion, visitedNeighbours);                
            }
        }        
    }
    
    void GetRadiiWithinDistance(unsigned vertexId, const SearchRegion& searchRegion, std::set<unsigned>& visitedNeighbours, std::vector<double>& radii)
    {
        // Check if we have visited this vertex in the recursion
        if (visitedNeighbours.find(vertexId) != visitedNeighbours.end())
        {
            return;
        }
        // Mark vertex as visited, if not
        visitedNeighbours.insert(vertexId);
        
        // Store the radius
        radii.push_back(vertexRadius[vertexId]);
        
        // Recursively visit neighbouring nodes contained in the search region
        vtkSmartPointer<vtkAdjacentVertexIterator> neighbour_iterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        networkGraph->GetAdjacentVertices(vertexId, neighbour_iterator);
        while (neighbour_iterator->HasNext())
        {
            unsigned neighIndex = neighbour_iterator->Next();
            
            if (searchRegion.IsPointContained(vertexCoordinates[neighIndex]))
            {
                GetRadiiWithinDistance(neighIndex, searchRegion, visitedNeighbours, radii);
            }
        }        
    }    
    
    vtkSmartPointer<vtkPolyData> networkSkeleton;
    vtkSmartPointer<vtkMutableUndirectedGraph> networkGraph;    
    std::vector<double*> vertexCoordinates;
    std::vector<double> vertexRadius;
    unsigned numVertices;
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    // Check for proper number of arguments
    if (nrhs!=3) 
    { 
	    mexErrMsgTxt("Three input arguments required."); 
    } 
    else if (nlhs!=0) 
    {
	    mexErrMsgTxt("No output arguments should be expected."); 
    } 

    // Get access to MATLAB data through C data types   
    std::string skeleton_filename(mxArrayToString(prhs[0]));
    std::string output_filename(mxArrayToString(prhs[1]));

    double* degree_3_crossings = mxGetPr(prhs[2]);
    unsigned num_degree_3_crossings = mxGetM(prhs[2]);
    if ((num_degree_3_crossings > 0) && (mxGetN(prhs[2])!=2))
    { 
	    mexErrMsgTxt("Third argument must be a Nx2 matrix with pairs of nodes defining pair-of-degree-3-type crossings"); 
    }

    NetworkAnalyser analyser(skeleton_filename);
    
    analyser.FindAndCorrectDegree4VesselCrossings();
    
    std::vector<std::pair<unsigned, unsigned> > crossing_list;
    
    for (unsigned crossing_id=0; crossing_id<num_degree_3_crossings; ++crossing_id)
    {
        crossing_list.push_back(std::pair<unsigned, unsigned>(degree_3_crossings[crossing_id], degree_3_crossings[num_degree_3_crossings+crossing_id])); 
    }
    
    analyser.Correct2PairsOfDegree3VesselCrossings(crossing_list);
    
    analyser.WriteToDisc(output_filename);
}
