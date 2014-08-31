/*
 *  dependencyGraph.h
 *  warpm Xcode project
 *
 *  Created by Noah Reddell on 5/19/11.
 *  Copyright 2011 University of Washington. All rights reserved.
 *
 */
#ifndef wmdependencygraph_h
#define wmdependencygraph_h

#include <utility> // std::pair
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/thread/mutex.hpp>

#include <iostream>
#include <fstream>

#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxstepperstatus.h>


/** \addtogroup lib 
 *  @{
 */

/**
 * Represents enumerated types of compute region completion and partial completion.
 *
 **/
class SubregionCompletionType
{
public:
	
	static const char* edgeTypeEnumText[];
	static unsigned numEdgeTypes;
	
	SubregionCompletionType()
	: completion_needed(edgeTypeEnumText[0]) //default to 'all' completion type
	{  
	};
	
	SubregionCompletionType( const std::string & desiredType)
	{  
		const char * matchedType;
		if (lookupCompletionType( desiredType, matchedType ) )
		{
			completion_needed = matchedType;
			return;
		}
		WxExcept wxe;
		wxe << "WmDependencyGraph<VertexProperty>::SubregionCompletionType() :Requested type not one of enumerated allowed types."  << std::endl;
		throw wxe;  // desired type does not match those possible
	};
	
	
	
	/**
	 * Sets the completion type property based on matching string specification
	 * @param desiredType should match the enumerated list of completion types possible
	 * @return true on success, false on failure (type was bad, no change made).
	 **/
	bool setCompletionType( const std::string & desiredType )
	{
		const char * matchedType;
		if (lookupCompletionType( desiredType, matchedType ) )
		{
			completion_needed = matchedType;
			return true;
		}
		else
			return false; // desired type does not match those possible
	};
	
	/**
	 * returns a string describing the type of completion needed
	 * 
	 * @return String associated with the completion type specified.
	 **/
	std::string getCompletionType() const
	{
		return std::string(completion_needed);
	};
	
	/**
	 * checks if right hand side completion type matches this one OR is a superset / more encompassing.
	 * as implemented now, the superset is simply the 'all' completion type, which trumps all others.
	 * 
	 * @return true if right hand side matches or is superset.
	 **/
	
	bool operator<=( const SubregionCompletionType & rhs ) const
	{
		if (rhs.completion_needed == edgeTypeEnumText[0] /* 'all' */)
			return true;
		if (rhs.completion_needed == completion_needed)
			return true;
		else
			return false;
		
	};
	
private:
	bool lookupCompletionType( const std::string & desiredType, const char * & found ) const
	{
		for (unsigned i = 0; i<numEdgeTypes; i++) {
			if (desiredType == edgeTypeEnumText[i] ) {
				found = edgeTypeEnumText[i];
				return true;
			}
		}
		return false;  // desired type does not match those possible
	};
	const char * completion_needed;
	
};


/**
 * Represents a graph relationship of tasks that can represent serial dependency and
 * opportunity or parallelism.
 *
 **/
template <class VertexProperty, class EdgeProperty>
class WmDependencyGraph
{
public:
	
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty, EdgeProperty > GraphType;
	typedef typename boost::graph_traits< GraphType >::vertex_descriptor VertexType;
	typedef typename boost::graph_traits< GraphType >::edge_descriptor   EdgeType;
	
	WmDependencyGraph()
	: graphStructureIsFixed(false)
	{
		
	};
	
	/**
	 * adds a new execution step vertex to the graph with specified prerequisite execution steps
	 * (zero is ok).
	 * @param newStep properties of the new step that will be copied
	 * @param preReqs list of prerequisite steps, and the type of prereq
	 * @return Vertex descriptor of new step is returned.
	 **/
	VertexType addExecutionStep( const VertexProperty & newStep,  const std::list<  std::pair< VertexType, EdgeProperty > > preReqs )
	{
		if (graphStructureIsFixed) {
			WxExcept wxe;
			wxe << "WmDependencyGraph<VertexProperty>::addExecutionStep : graph structure is fixed. No modifications allowed. ";
			throw wxe;
		}
		
		VertexType v = boost::add_vertex(boostgraph);
		
		// use default assignment operator to copy the vertex properties
		(boostgraph)[v] = newStep;
		
		bool b; EdgeType e;
		//insert into graph with with in-edges as specified by preReqs
		for(typename std::list<std::pair< VertexType, EdgeProperty > >::const_iterator list_iter = preReqs.begin(); 
		    list_iter != preReqs.end(); list_iter++)
		{
			tie(e, b) = add_edge((*list_iter).first, v, boostgraph);
			if (!b) {
				WxExcept wxe;
				wxe << "WmDependencyGraph<VertexProperty>::addExecutionStep : bad edge add ";
				throw wxe;
				return v;
			}
		}
		return v;
	};
	
	
	/**
	 * inserts a linear sequence graph into this one, making an append on the out-edge of the vertex specified
	 * verification is made that the supplied graph is linear.
	 * @param linear graph that will be copied and appended to this one
	 * @param vertex descriptor for the point at which graph will be appended (as an out edge)
	 * @return Vertex descriptor of tail of new sequence in the graph is returned.
	 **/
	VertexType appendSequence( const WmDependencyGraph<VertexProperty, EdgeProperty> & sourceSequence,  VertexType appendPoint, EdgeProperty edgetype )
	{
		if (graphStructureIsFixed) {
			WxExcept wxe;
			wxe << "WmDependencyGraph<VertexProperty>::appendSequence : graph structure is fixed. No modifications allowed. ";
			throw wxe;
		}
		
		VertexType tail, insertionPoint, sourceVertex;
		// first ensure that the supplied new sequence graph is linear
		if (!sourceSequence.isLinearGraph()) {
			WxExcept wxe;
			wxe << "WmDependencyGraph<VertexProperty>::addExecutionStep : source sequence is not linear";
			throw wxe;
			return tail;
		}
		
		typename boost::graph_traits<GraphType >::vertex_iterator sourcei, source_end;
		insertionPoint = appendPoint;
		
		//find head of new sequence, sourcei points to the head once out of loop
		for (boost::tie(sourcei, source_end) = boost::vertices(sourceSequence); sourcei != source_end; ++sourcei)
		{
			int id = in_degree( *sourcei, sourceSequence );
			if ( id == 0 )
				break;        
		}
		sourceVertex = *sourcei;
		
		// add each vertex in new sequence until reaching end.
		typename boost::graph_traits<GraphType>::out_edge_iterator sourceoei, sourceoe_end;
		
		int numSourceVertices = num_vertices( sourceSequence );
		for (int i= 0; i<numSourceVertices; i++) {
			std::list<VertexType> predecessor(1, insertionPoint );
			
			insertionPoint = this->addExecutionStep( sourceSequence[ sourceVertex],  predecessor);
			//advance to next vertex in source
			boost::tie(sourceoei, sourceoe_end) = boost::out_edges(sourceVertex, sourceSequence);
			
			if (i == numSourceVertices - 1) {
				break;  // below target call seg faults if no out edge.
			}
			sourceVertex = target( *sourceoei, sourceSequence );
		}
		
		
		
		tail = insertionPoint;
		return tail;
		
	};
	
	std::list< std::pair< VertexType, VertexProperty  >  >   getInDegreeZeroTasks() const
	{
		std::list< std::pair< VertexType, VertexProperty  >  >  returnList;
		typename boost::graph_traits<GraphType >::vertex_iterator vi, vi_end;
		for (boost::tie(vi, vi_end) = boost::vertices(boostgraph); vi != vi_end; ++vi)
		{
			if ( in_degree( *vi, boostgraph ) == 0 )
				returnList.push_back(  std::pair< VertexType, VertexProperty  > ( *vi, boostgraph[*vi]    )   );
		}
		return returnList;
	};
	
	std::list< std::pair< VertexType, VertexProperty  >  >   getAllVerticesAndTasks() const
	{
		std::list< std::pair< VertexType, VertexProperty  >  >  returnList;
		typename boost::graph_traits<GraphType >::vertex_iterator vi, vi_end;
		for (boost::tie(vi, vi_end) = boost::vertices(boostgraph); vi != vi_end; ++vi)
		{
			returnList.push_back(  std::pair< VertexType, VertexProperty  > ( *vi, boostgraph[*vi]    )   );
		}
		return returnList;
	};
	
	unsigned getNumWithOutDegreeZero() const
	{
		unsigned count = 0;
		typename boost::graph_traits<GraphType >::vertex_iterator vi, vi_end;
		for (boost::tie(vi, vi_end) = boost::vertices(boostgraph); vi != vi_end; ++vi)
		{
			if ( out_degree( *vi, boostgraph ) == 0 )
				count++;
		}
		return count;
	};
	
	std::list< VertexType > getDependentVertices( VertexType sourceVertex, EdgeProperty edgeprop) const
	{
		std::list< VertexType >  returnList;
		typename boost::graph_traits<GraphType>::out_edge_iterator sourceoei, sourceoe_end;
		for (boost::tie(sourceoei, sourceoe_end) = boost::out_edges(sourceVertex, boostgraph); sourceoei != sourceoe_end; sourceoei++ )
		{
			// check that the edge type matches
			if (boostgraph[*sourceoei] <= edgeprop ) {
				returnList.push_back( target(*sourceoei, boostgraph ) );
			}
		}
		
		return returnList;
		
	}
	
	const VertexProperty & getVertexProp( VertexType vertex ) const
	{
		return boostgraph[vertex];
	}
	
	unsigned getInDegree( VertexType vertex ) const
	{
		return in_degree( vertex, boostgraph );
	}
	
	
	/**
	 * a linear graph is one where all vertexes are serially chained together, so there is one vertex
	 * with in-degree zero, one with out-degree zero, and all others have in&out degree of one.
	 * @return true if linear graph criteria met.
	 **/
	bool isLinearGraph() const
	{
		/**
		 * one other case is possible that would meet the tested criteria:  One vertex with no in/out edges and all other
		 * nodes in a serial cyclic ring.  Since the WmDependencyGraph interface does not allow adding edges after vertex
		 * insertion, cyclic graphs are not possible, and thus this case does not need to be tested.
		 **/
		
		typename boost::graph_traits<GraphType >::vertex_iterator vi, vi_end;
		int numWithInDegreeZero = 0;
		int numWithOutDegreeZero = 0;
		for (boost::tie(vi, vi_end) = boost::vertices(boostgraph); vi != vi_end; ++vi)
		{
			int od = out_degree( *vi, boostgraph );
			if ( od == 0 )
				numWithOutDegreeZero++;
			else if ( od != 1 )
				return false;
			
			int id = in_degree( *vi, boostgraph );
			if ( id == 0 )
				numWithInDegreeZero++;
			else if ( id != 1 )
				return false;
			
		}
		
		if ( ( numWithInDegreeZero != 1) || ( numWithOutDegreeZero != 1) )
			return false;
		else
			return true;    
	};
	
	
	/**
	 * to better understand graph behavior, this will output the graph in graphviz format to the file specified.
	 **/
	void writeDebugOut(const std::string filename) const
	{
		std::ofstream out;
		out.open(filename.c_str() );
		write_graphviz(out, boostgraph ); //, make_label_writer(name));   
		out.close();
	};
	
	/**
	 * set to prevent subsequent modifications of this graph structure
	 **/
	void fixGraph()
	{
		graphStructureIsFixed = true;
	};
	
	/**
	 * true if graph is fixed to prevent subsequent modifications of this graph structure
	 **/
	bool isFixed() const
	{
		return graphStructureIsFixed;
	};
	
	
private:
	
	GraphType boostgraph;
	
	bool graphStructureIsFixed;
	
};

/** @}*/
#endif // wmdependencygraph_h