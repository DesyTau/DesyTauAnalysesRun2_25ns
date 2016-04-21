#ifndef KAPPA_RUNLUMIREADER_H
#define KAPPA_RUNLUMIREADER_H

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

//#include <Kappa/DataFormats/interface/Kappa.h>
//#include "../RootTools/FileInterface.h"

typedef unsigned int run_id;
typedef unsigned int lumi_id;

void readLumiFilter(const std::string json, std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > > &lumifilter);

class RunLumiSelector
{
public:
	RunLumiSelector(const std::string json = "", const run_id _passRunLow = 0, const run_id _passRunHigh = 0);
	RunLumiSelector(const std::vector<std::string> &json, const run_id _passRunLow = 0, const run_id _passRunHigh = 0);
	void addJSONFile(const std::string json = "");
	inline bool accept(const run_id run, const lumi_id lumi) const
	{
		if ((run >= passRunLow) && (passRunHigh > 0) && (run <= passRunHigh)) return(true);
		
		typedef std::set<std::pair<lumi_id, lumi_id> > lumirange;
		std::map<run_id, lumirange>::const_iterator itRun = lumifilter.find(run);
		if (itRun == lumifilter.end())
			return false;
		const lumirange &lumis = itRun->second;
		for (lumirange::const_iterator itLumis = lumis.begin(); itLumis != lumis.end(); ++itLumis)
			if (lumi < itLumis->first)
				return false;
			else
				if ((lumi >= itLumis->first) && (lumi <= itLumis->second))
					return true;
		return false;
	}
	std::pair<run_id, run_id> getBoundaries() const;
	std::pair<run_id, lumi_id> getMaxRunLumiPair() const;

	void printJSON(std::ostream &os = std::cout) const;
	static void printJSON(std::ostream &os,
		const std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > > &lumiSrc);
	static std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > >
		getMinimalJSON(const std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > > &lumiSrc);

	const std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > > &getRunLumiMap() const;
	friend std::ostream &operator<<(std::ostream &os, const RunLumiSelector &m);
private:
	run_id passRunLow, passRunHigh;
	std::map<run_id, std::set<std::pair<lumi_id, lumi_id> > > lumifilter;
};

std::ostream &operator<<(std::ostream &os, const RunLumiSelector &m);

#endif
