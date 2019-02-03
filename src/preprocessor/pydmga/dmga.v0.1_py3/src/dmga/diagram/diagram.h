/*
 * diagram.h
 *
 *  Created on: 11-11-2012
 *      Author: Robson
 */

#ifndef DIAGRAM_H_
#define DIAGRAM_H_

#include <vector>
#include <model/model.hpp>
#include <utils/utils.hpp>
#include <algorithm>
#include <algebra/algebra.hpp>
#include <algebra/polynomial.hpp>

#include <base.h>
#include <3rd/voro/voro++.hh>
#include <utils/voropp_container_relocator.hpp>


namespace dmga{

namespace diagram{

using namespace algebra;
using namespace polynomial;

/**
 * TODO: this class is supposed to be an interface for generic DIagram that computes Voronoi Diagram in 3D
 */
class GenericDiagram : public DMGAObject{
public:
	std::string show(){ return "dmga::diagram::GenericDiagram"; }

	~GenericDiagram(){ 			DEB3("GenericDiagram::__destruct__():");
		//destructor...
	}
};
/**
 * TODO: this class is supposed to be an interface for generic Kinetic Diagram that computes Voronoi Diagram in 3D and
 * allows for its update with time.
 */
class GenericKineticDiagram : public DMGAObject{
public:
	std::string show(){ return "dmga::diagram::GenericKineticDiagram"; }

	~GenericKineticDiagram(){ 			DEB3("GenericKineticDiagram::__destruct__():");
		//destructor...
	}
};

/**
 * Basic implementation of Voronoi Diagram using routines from Voro++
 * it gives nice wrapping for them as they are rather 'closed'.
 *
 * This class strongly depends on Voro++ implementation and uses its non-public api
 * (internal representation),
 */
template<typename ContainerSpec, typename CellSpec, typename StorageSpec = voro::container_relocator /*voro::container_poly*/>
class VoroppDiagram : public GenericDiagram {
public:
	std::string show(){ return "dmga::diagram::VoroppDiagram"; }
	typedef ContainerSpec ContainerType;
	typedef typename ContainerType::GeometryType GeometryType;
	typedef typename ContainerType::ElementType ParticleType;
	typedef typename ContainerType::KeyType KeyType;
	typedef CellSpec CellType;
	typedef std::vector<CellType*> ResultSetType;
	typedef std::vector<KeyType*> SubsetType;
	typedef std::vector<CellType*> CacheType;
	typedef StorageSpec StorageType; //TODO: narazie globalnie, byc moze lepiej zrobic tylko w KineticDiagram

	typedef int CachePolicyType;
	static const CachePolicyType CACHE_ON;
	static const CachePolicyType CACHE_OFF;

	/**
	 * this is location element, that holds
	 * information of location of particle inside
	 * voro++ container
	 */
	class LocationElement{
	public:
		int ijk;
		int q;
	};
	typedef std::vector<LocationElement> LocationType;
	/**
	 * this is Container class that can hold Particle information
	 * we use reference as we do not allow to create diagram without container!
	 */
	ContainerType& container;
	/**
	 * here is the voro++ container to do computations
	 * we need pointer because pointer is something that
	 * ContainerType produces
	 */
	StorageType* storage;
	/**
	 * here we will store corespondence between cell number ad its location inside
	 * the voro++ container (it's this (ijk, q) pair)
	 * we need pointer because pointer is something that
	 * ContainerType produces
	 */
	LocationType* location;
	/**
	 * this is cache that holds computed cells if CACHE_ON is used
	 */
	CacheType cache;
	/**
	 * this is workaround for the problem of too big voronoi diagrams!
	 * if caching is off then no Cell is stored inside Diagram, only
	 * recomputed every time, the user is responsible for freeing them
	 *
	 * if the cache is ON then Diagram will free all alocated memory upon destruction
	 *
	 * TODO: przemyslec?
	 * */
	CachePolicyType cache_policy;


	/**
	 * make new Diagram for container
	 * if cache is on then it compytes all cells upon creation
	 */
	VoroppDiagram(ContainerType& container, CachePolicyType cache_policy = CACHE_ON) :
			container(container),
			storage(0),
			location(0),
			cache(container.size()),
			cache_policy(cache_policy){																						DEB3("VoroppDiagram::__construct__()");
		int n = container.size();
		container.makeVoroppStorage(storage, location);
		if (cache_policy == CACHE_ON){
			for (int i = 0; i < n; ++i){ 																					DEB3("VoroppDiagram::__construct__() computing " << i << "-th cell...");
				cache[i] = new CellType(i); 																				DEB3("VoroppDiagram::__construct__() cell created");
				computeCell(i);																								DEB3("VoroppDiagram::__construct__() cell computed");
			}
		}
	}

	/**
	 * this constructor assumes CACHE_POLICY_ON
	 *
	 * upon creation only a subset of Voronoi diagram is created
	 */
	VoroppDiagram(ContainerType& container, const SubsetType& subset) :
			container(container),
			storage(0),
			location(0),
			cache(container.size()),
			cache_policy(CACHE_ON){																							DEB3("VoroppDiagram::__construct__(): with subset...");
		container.makeVoroppStorage(storage, location);
		auto it = subset.begin();
		auto end = subset.begin();
		for (; it < end; ++it){
			cache[(int)(*it)] = new CellType((int)(*it));
			computeCell((int)(*it));
		}
	}
	/**
	 * returns i-th Cell as a pointer
	 */
	CellType* getCell(int i){																								DEB3("VoroppDiagram::getCell(): getting " << i);
		return getCachedCell(i);
	}
	/**
	 * return i-th Cell as a refernce
	 *
	 * TODO: dont use it when CACHE_POLICY_OFF! Lost pointer!
	 */
	CellType& operator[](int i){																							DEB3("VoroppDiagram::operator[](): getting " << i);
		return *getCachedCell(i);
	}
	/**
	 * return all Cells
	 *
	 * if necessary it computes missing Cells
	 */
	ResultSetType getCells(){																								DEB3("VoroppDiagram::getCells()");
		int n = container.size();
		ResultSetType result(n);
		for (int i = 0; i < n; ++i){
			result[i] = getCachedCell(i);
		}
		return result;
	}
	/**
	 * returns only a subset of cells
	 *
	 * if necessary it computes missing Cells
	 */
	ResultSetType getCells(const SubsetType& subset){																		DEB3("VoroppDiagram::getCells(): with subset...");
		throw dmga::exceptions::NotImplementedYet("VoroppDiagram::getCells()");
	}
	/**
	 * destructor - frees all Cells in cache
	 */
	~VoroppDiagram(){																										DEB3("VoroppDiagram::__destruct__");
		for (auto it = cache.begin(); it != cache.end(); ++it){
			dmga::utils::safeDelete(*it);
		}
		dmga::utils::safeDelete(storage);
		dmga::utils::safeDelete(location);
	}
	/**
	 * This is experimental
	 *
	 * Adds preset to this container that allows for initialization
	 * of a cell before it is computed, that may be usefull for
	 * nonreagular packaging
	 */
	void addPreset(unsigned int i, dmga::model::BaseCellPreset& preset){
		assureCachedCell(i);
		preset.setTarget(&(this->cache[i]->voropp_cell));
		DEB2("pre add preset");
		storage->add_wall(preset);
		DEB2("post add preset");
	}
protected:
	/**
	 * recomputes i-th cell
	 */
	CellType* computeCell(unsigned int i){																					DEB3("VoroppDiagram::computeCell(): computing " << i);
		int ijk = (*location)[i].ijk;
		int q = (*location)[i].q;
		if(storage->compute_cell(cache[i]->voropp_cell, ijk, q)){															DEB3("VoroppDiagram::computeCell(): voro++ computed");
			cache[i]->setEmpty(false);
			voro::normalizeCell(cache[i]->voropp_cell, container.get(i).getRaw()/*dmga::utils::getContainerCoordsRaw(*storage, ijk, q)*/);				DEB3("VoroppDiagram::computeCell(): normalized");
			cache[i]->update(); 																							DEB3("VoroppDiagram::computeCell(): updated");
		}else{
			cache[i]->setEmpty(true);
		}
		return cache[i];
	}
	/**
	 * assures that i-th cell has place inside cache container
	 */
	void assureCachedCell(unsigned int i){																					DEB3("VoroppDiagram::resizeCache(): new size " << i);
		if (cache.size() <= i){
			while (cache.size() < i){
				cache.push_back(0);
			}
			cache.push_back(new CellType(i));
		}else if (!cache[i]){
			cache[i] = new CellType(i);
		}
	}
	/**
	 * return cached cell
	 *
	 * if CACHE_POLICY_OF then it allways recomputes cell from scratch
	 */
	CellType* getCachedCell(unsigned int i){																				DEB3("VoroppDiagram::getCachedCell(): getting " << i);
		if (cache_policy == CACHE_ON){
			if ( !((cache.size() > i) && cache[i]) ){
				assureCachedCell(i);
				return computeCell(i);
			}
			return cache[i];
		}else{
			CellType* cell = new CellType(i);
			int ijk = (*location)[i].ijk;
			int q = (*location)[i].q;
			if (storage->compute_cell(cell->voropp_cell, ijk, q)){											 				DEB3("VoroppDiagram::computeCell(): voro++ computed");
				cell->setEmpty(false);
				voro::normalizeCell(cell->voropp_cell, container.get(i).getRaw()/*dmga::utils::getContainerCoordsRaw(*storage, ijk, q)*/);				DEB3("VoroppDiagram::computeCell(): normalized");
				cell->update(); 																							DEB3("VoroppDiagram::computeCell(): updated");
			}else{
				cell->setEmpty(true);
			}
			return cell;
		}
	}
};

template<typename ContainerSpec, typename CellSpec, typename StorageSpec>
const typename VoroppDiagram<ContainerSpec, CellSpec, StorageSpec>::CachePolicyType VoroppDiagram<ContainerSpec, CellSpec, StorageSpec>::CACHE_ON = 1;

template<typename ContainerSpec, typename CellSpec, typename StorageSpec>
const typename VoroppDiagram<ContainerSpec, CellSpec, StorageSpec>::CachePolicyType VoroppDiagram<ContainerSpec, CellSpec, StorageSpec>::CACHE_OFF = 0;

}//namespace diagram

}//namespace dmga


#endif /* DIAGRAM_H_ */
