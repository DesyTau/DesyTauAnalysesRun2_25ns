#ifndef DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_CONFIG_H
#define DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_CONFIG_H

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include "TDirectory.h"
#include "TObjString.h"

inline std::string trim(const std::string& str)
{
	std::string::size_type first = str.find_first_not_of(" \t\n\v");
	std::string::size_type last = str.find_last_not_of(" \t\n\v");

	if(first == std::string::npos || last == std::string::npos)
		return "";

	return str.substr(first, last - first + 1);
}

namespace ConfigConv
{
	template<typename T>
	struct conv
	{
		T operator()(const std::string& value) const
		{
			std::stringstream stream(value);
			T val;
			stream >> val;
			if(!stream || stream.bad())
				throw std::runtime_error("Failed to parse value " + value);
			return val;
		}
	};

	template<>
	struct conv<bool>
	{
		bool operator()(const std::string& value) const
		{
			if(value == "true" || value == "True") return true;
			if(value == "false" || value == "False") return false;
			throw std::runtime_error("Failed to parse boolean value " + value);
		}
	};

	template<>
	struct conv<std::string>
	{
		const std::string& operator()(const std::string& value) const { return value; }
	};

	template<typename T>
	struct conv<std::vector<T> >
	{
		std::vector<T> operator()(const std::string& value) const
		{
			if(value.empty()) return std::vector<T>();

			std::vector<T> result;
			std::string::size_type prev = 0;
			while(prev != std::string::npos)
			{
				const std::string::size_type pos = value.find(',', prev);
				if(pos == std::string::npos)
				{
					result.push_back(conv<T>()(trim(value.substr(prev))));
					prev = pos;
				}
				else
				{
					result.push_back(conv<T>()(trim(value.substr(prev, pos - prev))));
					prev = pos + 1; // skip ','
				}
			}

			return result;
		}
	};
}

class Config
{
public:
	Config(const char* filename);

	void merge(const Config& other);

	template<typename T>
	T get(const std::string& entry) const
	{
		std::map<std::string, std::string>::const_iterator iter = entries.find(entry);
		if(iter == entries.end()) throw std::runtime_error("Failed to find configuration entry " + entry);

		return ConfigConv::conv<T>()(iter->second);
	}

	template<typename T>
	T get(const std::string& entry, T stdVal) const
	{
		std::map<std::string, std::string>::const_iterator iter = entries.find(entry);
		if(iter == entries.end()) return stdVal;

		return ConfigConv::conv<T>()(iter->second);
	}

	void writeConfigToTree()
	{
		for (std::map<std::string, std::string>::const_iterator it = entries.begin(); it != entries.end(); ++it)
		{
			TObjString tmpString = TObjString(it->second.c_str());
			tmpString.Write(it->first.c_str());
		}
	}
private:
	std::map<std::string, std::string> entries;
};

#endif // DESYHTAUTAU_TAUTAUANALYSIS_INTERFACE_CONFIG_H
