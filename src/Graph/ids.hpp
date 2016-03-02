/*
 * ids.hpp
 *
 *  Created on: Mar 1, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRAPH_IDS_HPP_
#define SRC_GRAPH_IDS_HPP_

/*! Here we define different the remapped-id
 *
 * rid, gid and lid are all unsigned long integer, and can be easily interchanged by mistake
 *  encapsulating avoid that this could happen. The second is readability, from the definition
 *  of function/structure we see immediately which id parameter accept/store
 *
 */
struct rid
{
	idx_t id;

	inline bool operator<=(const rid & r) const
	{
		return id <= r.id;
	}

	inline bool operator<(const rid & r) const
	{
		return id < r.id;
	}

	inline rid operator-(int i) const
	{
		struct rid tmp;
		tmp.id = id - i;
		return tmp;
	}

	inline rid operator-(struct rid i) const
	{
		struct rid tmp;
		tmp.id = id - i.id;
		return tmp;
	}

	inline rid operator+(int i) const
	{
		struct rid tmp;
		tmp.id = id + i;
		return tmp;
	}

	inline rid & operator+=(const rid & i)
	{
		id += i.id;
		return *this;
	}

	inline rid & operator++()
	{
		id++;

		return *this;
	}

	inline bool operator==(const rid & r) const
	{
		return id == r.id;
	}
};

/*! Here we define different the remapped-id
 *
 * rid, gid and lid are all unsigned long integer, and can be easily interchanged by mistake
 *  encapsulating avoid that this could happen. The second is readability, from the definition
 *  of function/structure we see immediately which id parameter accept/store
 *
 */
struct gid
{
	size_t id;
};

/*! Here we define different the remapped-id
 *
 * rid, gid and lid are all unsigned long integer, and can be easily interchanged by mistake
 *  encapsulating avoid that this could happen. The second is readability, from the definition
 *  of function/structure we see immediately which id parameter accept/store
 *
 */
struct lid
{
	size_t id;
};

// define hash map for gid rid and lid

namespace std
{
	template <>
	struct hash<rid>
	{
		inline std::size_t operator()(const rid& k) const
		{
			return k.id;
		}
	};

	template <>
	struct hash<gid>
	{
		inline std::size_t operator()(const gid& k) const
		{
			return k.id;
		}
	};

	template <>
	struct hash<lid>
	{
		inline std::size_t operator()(const lid& k) const
		{
			return k.id;
		}
	};

}


#endif /* SRC_GRAPH_IDS_HPP_ */
