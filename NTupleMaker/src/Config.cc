#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"

#include <fstream>
#include <cstring>

Config::Config(const char* filename)
{
	std::ifstream stream(filename);
	if(!stream)
		throw std::runtime_error(std::string("Failed to open configuration file ") + filename);

	std::string line;
	while(std::getline(stream, line))
	{
		// filter comments and empty lines
		std::string::size_type com_pos = line.find_first_of("#;");
		if(com_pos != std::string::npos) line.erase(com_pos);
		line = trim(line);
		if(line.empty()) continue;

		std::string::size_type pos = line.find('=');
		if(pos == std::string::npos)
			throw std::runtime_error("Invalid line in configuration file: " + line);

		std::string key = trim(line.substr(0, pos));
		std::string value = trim(line.substr(pos+1));

		if(key == "include")
		{
			const char* pos = std::strrchr(filename, '/');
			if(pos != NULL)
				value = std::string(filename, pos - filename + 1) + value;
			const Config included(value.c_str());
			for(std::map<std::string, std::string>::const_iterator iter = included.entries.begin(); iter != included.entries.end(); ++iter)
				if(entries.count(iter->first) == 0) entries[iter->first] = iter->second;
		}
		else
		{
			entries[key] = value;
		}
	}
}

void Config::merge(const Config& other)
{
	for(std::map<std::string, std::string>::const_iterator iter = other.entries.begin(); iter != other.entries.end(); ++iter)
	{
		std::map<std::string, std::string>::const_iterator my_iter = entries.find(iter->first);
		if(my_iter != entries.end() && iter->second != my_iter->second)
			throw std::runtime_error("Merge: Different values for key=" + iter->first);
		entries[iter->first] = iter->second;
	}
}
