/*
 * container.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <base.h>
#include <exceptions/exceptions.hpp>
#include <3rd/voro/voro++.hh>


namespace dmga{

namespace container{

/**
 * defines generic type of container, givin an interface, that
 * any container should have...
 */
template<typename GeometrySpec, typename ElementSpec>
class GenericContainer : public DMGAObject{
public:
	std::string show(){ return "dmga::container::GenericContainer"; }
	typedef GeometrySpec GeometryType;
	typedef ElementSpec ElementType;
	typedef typename ElementType::KeyType KeyType;

	GeometryType& geometry;

	/**
	 * dodaje element do pojemnika
	 */
	virtual void add(const ElementType& element) = 0;
	/**
	 * zwraca i-ty wstawiony element
	 */
	virtual ElementType get(int i) = 0;
	/**
	 * znajduje element o podanym kluczu
	 */
	virtual ElementType find(const KeyType& key) = 0;

	virtual ElementType& operator[](int index){
		throw dmga::exceptions::NotSupported("This container doesn't support mutable elements");
	}
	/**
	 * tworzy nowy pojemnik - geometria jest wymagana!
	 */
	GenericContainer(GeometryType& geometry) : geometry(geometry){		DEB3("GenericContainer::__construct__(): ");
	}

	virtual ~GenericContainer(){										DEB3("GenericContainer::__destruct__(): ");
	}
};

/**
 * Generic Voro++ compatibile Container -
 * those functions sould be overriden to use standard Voro++ Diagram types
 */
template<typename GeometrySpec, typename ElementSpec >
class GenericVoroppContainer : public GenericContainer<GeometrySpec, ElementSpec>{
public:
	typedef GenericContainer<GeometrySpec, ElementSpec> super;
	typedef GeometrySpec GeometryType;
	typedef ElementSpec ElementType;
	typedef typename ElementType::KeyType KeyType;

	GenericVoroppContainer(GeometryType& geometry) : super(geometry){						DEB3("GenericVoroppContainer::__construct__(): ");
	}

	template<typename StorageType, typename LocationType>
	void makeVoroppStorage(StorageType* some_storage, LocationType* some_location){			DEB3("GenericVoroppContainer::makeVoroppStorage(): ");
	}

	template<typename StorageType, typename LocationType>
	void updateVoroppStorage(StorageType* some_storage, LocationType* some_location){		DEB3("GenericVoroppContainer::updateVoroppStorage(): ");
	}

	virtual ~GenericVoroppContainer(){														DEB3("GenericVoroppContainer::__destruct__(): ");
	}
};

/**
 * Default simple Voro++ compatibile Container -
 */
template<typename GeometrySpec, typename ElementSpec>
class VectorContainer : public GenericVoroppContainer<GeometrySpec, ElementSpec> {
public:
	typedef GenericVoroppContainer<GeometrySpec, ElementSpec> super;
	typedef GeometrySpec GeometryType;
	typedef ElementSpec ElementType;
	typedef typename ElementType::KeyType KeyType;

	/** make new storage */
	template<typename StorageType, typename LocationType>
	void makeVoroppStorage(StorageType*& some_storage, LocationType*& some_location){		DEB3("VectorContainer::makeVoroppStorage(): ");
		//find the best partition of the box... (only very simple algorithm now... )
		int num_of_boxes = (size() + 6) / 6; //we assume that this is the best number of particles per box, we add +6 because then we always get num_of_boxes > 0
		int est = 1;
		//find first bigger
		while (est * est * est < num_of_boxes) est <<= 1; //est *= 2
		//two steps of newton operator...
		double d_est = est;
		d_est = d_est - ((d_est - ((double)num_of_boxes / (d_est * d_est))) / 3.0);
		est = (int)(d_est - ((d_est - ((double)num_of_boxes / (d_est * d_est))) / 3.0));
		//if 0 then 1
		if (est <= 0) est = 1;
		//TODO: przemyslec czy ta implementacja jest dobra :)

		this->geometry.makeBox(some_storage, est, est, est);
		some_location = new LocationType(data.size());
		auto data_it = data.begin();
		auto data_end = data.end();
		int num = 0;
		for (; data_it < data_end; ++data_it, ++num){
			dmga::utils::putContainerCoordsRaw(*some_storage, coords_fetch, num, (*data_it).getRaw());
			(*some_location)[num].ijk = coords_fetch.ijk;
			(*some_location)[num].q = coords_fetch.q;
		}
	}
	/** updates a storage */
	template<typename StorageType, typename LocationType>
	void updateVoroppStorage(StorageType*& some_storage, LocationType*& some_location){			DEB3("VectorContainer::updateVoroppStorage(): ");
		//this is old way, without container_relocator - we destroy container and create new one
		dmga::utils::safeDelete(some_storage);
		dmga::utils::safeDelete(some_location);
		makeVoroppStorage(some_storage, some_location);

		//this is a new way - using relocator - we suppose that this should be somehow faster
//		some_storage->batchset_start();
//		for (unsigned int i = 0; i < data.size(); ++i){
//			double* raw = data[i].getRaw();
//			some_storage->batchset_enqueue((*some_location)[i].ijk, (*some_location)[i].q, raw[0], raw[1], raw[2]);
//		}
//		some_storage->batchset_finalize(*some_location);
	}
	/** safely deletes storage */
	~VectorContainer(){																									DEB3("VectorContainer::__destruct__(): ");
	}
	/** creates new Container for given geometry */
	VectorContainer(GeometryType& geometry) : super(geometry){															DEB3("VectorContainer::__construct__(): ");											DEB3("VoroppContainer::__construct__()");
	}
	/** adds new particle */
	void add(const ElementType& element){																				DEB3("VectorContainer::add(): adding " << size() << " " << element.getX() << " " << element.getY() << " " << element.getZ() << " " << element.getR());
		data.push_back(element);
	}
	/** returns i-th particle (the one that was added as an i-th particle, counted from 0) */
	ElementType get(int i){																								DEB3("VectorContainer::get(): getting " << i);
		return data[i];
	}
	/** finds a particle for a given key type */
	ElementType find(const KeyType& key){																				DEB3("VectorContainer::find(): searching " << key);
		throw exceptions::NotImplementedYet("VectorContainer::find()");
	}
	/** returns current coords in raw format (array of doubles) */
	double* getRawCoords(int i){																						DEB3("VectorContainer::get(): getting " << i);
		return data[i].getRaw();
	}
	/** returns mutable element, you may need to update any diagram using this container */
	ElementType& operator[](int index){
		return data[index];
	}
	/** updates all particles coordinates using Trajectory class */
	template<typename TrajectorySpec>
	void update(TrajectorySpec& trajectory){																			DEB3("VectorContainer::update(): by some trajectory... ");
		auto data_it = data.begin();
		auto data_end = data.end();
		double* raw = trajectory.getRaw();
		for (; data_it < data_end; ++data_it){
			double* ptr = (*data_it).getRaw();
			double* last = (ptr + 3); //only 3 first coordinates, don't touch radius yet...
			while (ptr < last){
				*ptr = *raw;
				++ptr; ++raw;
			}
		}
	}
	/**
	 * returns current number of particles
	 */
	inline int size() const{																							DEB3("VectorContainer::size(): returns " << data.size());
		return data.size();
	}

protected:
	std::vector<ElementType> data;

	/** this class allows to remember the order of puting particles, with it we are able to retrive ijk and q indices for each particle */
	typedef voro::fetching_particle_order InternalCoordFetcherType;
	/** this will allow us to retrive ijk and q indices of the underlying voro++ implementation */
	InternalCoordFetcherType coords_fetch;
};

/**
 * Standard container for Voro++ library part
 * it uses voro::container_poly as a storage type
 * to hold particle positions
 *
 * ElementSpec tells what kind of elements we want to store in container,
 * it should at last provide 4 coordinates (x,y,z,r) and an id of some KeyType
 * For example dmga::model::WeightedParticle is a good example
 *
 * GeometrySpec tells which geometry class should be used to define box for
 * particles. For example dmga::geometry::OrthogonalGeometry
 *
 * this container type allows updating of particle positions by Trajectory
 *
 * @deprecated - change in the structure, use other Containers (as VectorContainer) as they
 *               are more object-oriented now and more generic
 */
template<typename GeometrySpec, typename ElementSpec>
class VoroppContainer : public GenericVoroppContainer<GeometrySpec, ElementSpec>{
public:
	typedef GenericVoroppContainer<GeometrySpec, ElementSpec> super;
	typedef typename super::GeometryType GeometryType;
	typedef typename super::ElementType ElementType;
	typedef typename super::KeyType KeyType;
	typedef typename voro::container_relocator StorageType; //specyfic voro++ container storage..., should be poly

private:
	/**
	 * class to store more information on Particles
	 */
	class ElementInternalData{
	public:
		KeyType key;
		int ijk;
		int q;
		/**
		 * creates new element type pointing to appropriate
		 * element in voro::container_poly (by ijk and q indices)
		 */
		ElementInternalData(KeyType key, int ijk, int q) :
				key(key),
				ijk(ijk),
				q(q){
		}
		/** return actual coordinates of this particle */
		double* getActualCoords(StorageType& storage){
			return (storage.p[ijk] + ElementType::RAW_SIZE * q);
		}
	};
	/** this class allows to remember the order of puting particles, with it we are able to retrive ijk and q indices for each particle */
	typedef voro::fetching_particle_order InternalCoordFetcherType;
	/** here we store all the additional data on particles (other than position, radius and number) */
	typedef std::vector<ElementInternalData> LocationType;
	/** this will allow us to retrive ijk and q indices of the underlying voro++ implementation */
	InternalCoordFetcherType coords_fetch;
public:
	/** here we store all particle coordinates */
	StorageType* storage;
	/** here we can search for key to locate particle */
	LocationType location;
	/** makes new storage and returns it, user is responsible for destroing it later */
	StorageType* makeStorage(GeometryType& geometry){																	DEB3("VoroppContainer::makeStorage(): ");
		StorageType* storage;
		geometry.makeBox(storage, 1, 1, 1);
		return storage;
	}
	/** make new storage */
	template<typename LocationType>
	void makeVoroppStorage(StorageType*& some_storage, LocationType*& some_location){
		some_storage = makeStorage(this->geometry);
		some_location = new LocationType(location.size());
		auto loc_it = location.begin();
		auto loc_end = location.end();
		int ijk, q;
		int id = 0;
		for (; loc_it < loc_end; ++loc_it, ++id){
			ijk = (*loc_it).ijk;
			q = (*loc_it).q;
			(*some_location)[id].ijk = ijk;
			(*some_location)[id].q = q;
			double* curr = dmga::utils::getContainerCoordsRaw(*storage, ijk, q);
			dmga::utils::putContainerCoordsRaw(*some_storage, id, curr);
		}
	}
	/** updates a storage */
	template<typename LocationType>
	void updateVoroppStorage(StorageType*& some_storage, LocationType*& some_location){
		//zakladamy ze storage istnieje i ze maja takie same wszystkie pola
		auto loc_it = location.begin();
		auto loc_end = location.end();
		int ijk, q;
		int id = 0;
		for (; loc_it < loc_end; ++loc_it, ++id){
			ijk = (*loc_it).ijk;
			q = (*loc_it).q;
			(*some_location)[id].ijk = ijk;
			(*some_location)[id].q = q;
			double* curr = storage->p[ijk] + ElementType::RAW_SIZE * q;
			double* other = some_storage->p[ijk] + ElementType::RAW_SIZE * q;
			for (int z = 0; z < ElementType::RAW_SIZE; ++z) *other++ = *curr++;
		}
	}
	/** safely deletes storage */
	~VoroppContainer(){																									DEB3("VoroppContainer::__destruct__(): ");
		dmga::utils::safeDelete(storage);
	}
	/** creates new Container for given geometry */
	VoroppContainer(GeometryType& geometry) : super(geometry){															DEB3("GenericContainer::__construct__(): ");											DEB3("VoroppContainer::__construct__()");
		this->storage = makeStorage(geometry);
	}
	/** adds new particle */
	void add(const ElementType& element){																				DEB3("VoroppContainer::add(): adding " << size() << " " << element.getX() << " " << element.getY() << " " << element.getZ() << " " << element.getR());
		//TODO: IMPORTANT: !!!!!!!!!!!!!!!!!!!11111 Sprzawdzać, czy sie dodało!!!!!!!!!!
		dmga::utils::putContainerCoordsRaw(*storage, coords_fetch, size(), element.getRaw());
		//storage->put(coords_fetch, size(), element.getX(), element.getY(), element.getZ(), element.getR());
		location.push_back(ElementInternalData(element.getKey(), coords_fetch.ijk, coords_fetch.q));
	}
	/** returns i-th particle (the one that was added as an i-th particle, counted from 0) */
	ElementType get(int i){																								DEB3("VoroppContainer::get(): getting " << i);
		return ElementType(location[i].key, location[i].getActualCoords(*storage));
	}
	/** finds a particle for a given key type */
	ElementType find(const KeyType& key){																				DEB3("VoroppContainer::find(): searching " << key);
		throw exceptions::NotImplementedYet("VoroppContainer::find()");
	}
	/** returns current coords in raw format (array of doubles) */
	double* getRawCoords(int i){																						DEB3("VoroppContainer::get(): getting " << i);
		return location[i].getActualCoords(*storage);
	}
	/** updates all particles coordinates using Trajectory class */
	template<typename TrajectorySpec>
	void update(TrajectorySpec& trajectory){																			DEB3("VoroppContainer::update(): by some trajectory... ");
		auto loc_it = location.begin();
		auto loc_end = location.end();
		int ijk, q;
		double* raw = trajectory.getRaw();
		for (; loc_it < loc_end; ++loc_it){
			ijk = (*loc_it).ijk;
			q = (*loc_it).q;
			double* ptr = storage->p[ijk] + ElementType::RAW_SIZE * q;
			double* last = (ptr + 3); //only 3 first coordinates, don't touch radius yet...
			while (ptr < last){
				*ptr = *raw;
				++ptr; ++raw;
			}
		}
	}
	/**
	 * returns current number of particles
	 */
	inline int size() const{																							DEB3("VoroppContainer::size(): returns " << location.size());
		return location.size();
	}
};

}//namespace container

}//namespace dmga


#endif /* CONTAINER_H_ */
