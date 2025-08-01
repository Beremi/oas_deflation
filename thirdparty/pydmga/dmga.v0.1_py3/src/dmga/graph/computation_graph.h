/*
 * computation_graph.h
 *
 *  Created on: 26-03-2013
 *      Author: Robson
 */

#ifndef COMPUTATION_GRAPH_H_
#define COMPUTATION_GRAPH_H_

#include <vector>
#include <deque>

namespace dmga{

namespace graph{

/**
 * this class allows to control sequence of computation of cells
 * i.e. it allows to compute cells in some given order, skip computation of some cells,
 * add new cells to compute queue as one computes other cells etc.
 *
 * It allows for example to compute only those cells that are needed for specyfic task
 * for example only external water molecules when one wants to find the surface of
 * the lipid bilayer, or similar tasks.
 */
template<typename DiagramSpec, typename QueueSpec, typename ExpanderSpec>
class ComputationGraph{
public:
	typedef DiagramSpec DiagramType;
	typedef QueueSpec QueueType;
	typedef ExpanderSpec ExpanderType;
	typedef typename DiagramType::CellType CellType;

	DiagramType& diagram;
	QueueSpec queue;
	ExpanderType expander;
	std::vector<int> flag;

	static const int FLAG_WHITE = 0;
	static const int FLAG_GRAY = 1;
	static const int FLAG_BLACK = 2;

	ComputationGraph(DiagramType& diagram) : diagram(diagram), flag(diagram.container.size(), FLAG_WHITE){
	}

	ComputationGraph(DiagramType& diagram, const ExpanderType& expander) :
			diagram(diagram),
			expander(expander),
			flag(diagram.container.size(), FLAG_WHITE){
	}

	void setExpander(const ExpanderType& new_expander){
		expander = new_expander;
	}

	inline void append(int cell_number){
		if (flag[cell_number] == FLAG_WHITE){
			queue.push(cell_number);
			flag[cell_number] = FLAG_GRAY;
		}
	}

	template<typename IteratorSpec>
	void append(const IteratorSpec& begin, const IteratorSpec& end){
		IteratorSpec it = begin;
		for (; it != end; ++it){
			append(*it);
		}
	}

	bool hasNext(){
		return !queue.isEmpty();
	}

	CellType* next(){
		if (queue.isEmpty()){
			return 0;
		}
		int cell_number = queue.pop();
		flag[cell_number] = FLAG_BLACK;
		CellType* cell = diagram.getCell(cell_number);
		if (!cell->isEmpty()){ //TODO: think! moze dac mozliwosc?
			std::vector<int> expand_to_cells = expander(cell);
			append(expand_to_cells.begin(), expand_to_cells.end());
		}
		return cell;
	}

};

/**
 * simple First In First Out queue that can be accepted by ComputationGraph
 */
class FIFOQueue{
public:
	std::deque<int> queue;

	inline bool isEmpty(){
		return queue.size() == 0;
	}

	inline void push(int smth){
		queue.push_back(smth);
	}

	inline int top(){
		return queue.front();
	}

	inline int pop(){
		int elem = queue.front();
		queue.pop_front();
		return elem;
	}
};

/**
 * simple Last In First Out (Stack) queue that can be accepted by ComputationGraph
 */
class LIFOQueue{
public:
	std::vector<int> queue;

	inline bool isEmpty(){
		return queue.size() == 0;
	}

	inline void push(int smth){
		queue.push_back(smth);
	}

	inline int top(){
		return queue.back();
	}

	inline int pop(){
		int elem = queue.back();
		queue.pop_back();
		return elem;
	}
};

/**
 * accepts all cells
 */
class AllExpander{
public:
	std::vector<int> to_expand;

	template<typename CellSpec>
	std::vector<int>& operator()(CellSpec* cell){
		to_expand.clear();
		voro::CellIterator it = cell->getCellIterator();
		while (!it.isFinished()){
			voro::EdgeDesc edge = it.current();
			to_expand.push_back(cell->getNeighbourId(edge.v, edge.j));
			while(!it.isMarked()){
				it.mark();
				it.forward();
			}
			it.jump();
		}
		return to_expand;
	}
};

/**
 * base class for simple accept expander
 *
 * it needs to override "accept" function for specyfic atom
 * to decide if the computation should proceed to this cell
 */
class AcceptExpander{
public:
	std::vector<int> to_expand;

	template<typename CellSpec>
	std::vector<int>& operator()(CellSpec* cell){
		to_expand.clear();
		voro::CellIterator it = cell->getCellIterator();
		while (!it.isFinished()){
			voro::EdgeDesc edge = it.current();
			int neighbourNumber = cell->getNeighbourId(edge.v, edge.j);
			if (this->accept(neighbourNumber)){
				to_expand.push_back(neighbourNumber);
			}
			while(!it.isMarked()){
				it.mark();
				it.forward();
			}
			it.jump();
		}
		return to_expand;
	}

	virtual bool accept(int cell_number) = 0;

	virtual ~AcceptExpander(){
	}
};

}//namespace graph

}//namespace dmga


#endif /* COMPUTATION_GRAPH_H_ */
