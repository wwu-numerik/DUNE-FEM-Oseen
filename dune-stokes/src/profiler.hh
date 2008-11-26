#ifndef PROFILER_HH_INCLUDED
#define PROFILER_HH_INCLUDED

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <ctime>

struct TimingData
{
	clock_t start;
	clock_t end;
	std::string name;
	TimingData( const std::string _name, const clock_t _start ):start( _start ),end( 0.0 ),name( _name ) {};
	TimingData():start( 0.0 ),end( 0.0 ),name( "blank" ) {};

};

typedef std::map<std::string,TimingData> DataMap;
typedef std::vector<DataMap> MapVector;

class Profiler
{

	public:
		Profiler() {}
		~Profiler() {};
		inline void StartTiming( const std::string funcname )
		{
			TimingData td( funcname,clock() );
			(m_timings[m_cur_run_num])[funcname] = td;
		}

		inline void StopTiming( const std::string funcname )
		{
			(m_timings[m_cur_run_num])[funcname].end = clock();
		}

        template < class CollectiveCommunication >
		void Output( CollectiveCommunication& comm, const int refineLevel, const long numDofs );

		inline void Reset( const int numRuns )
		{
			m_timings.clear();
			//m_timings.reserve( numRuns*2 );
			m_total_runs = numRuns;
			m_cur_run_num = 0;
			DataMap td;
            m_timings.push_back( td );
            m_l2_error = 0;
		}

		inline void AddCount( const int num)
		{
		    m_count[num] +=1;
		}

		inline void NextRun(const double error )
		{
            m_cur_run_num +=  m_cur_run_num < m_total_runs - 1 ?  1 : 0 ;
            DataMap td;
            m_timings.push_back( td );
            m_l2_error += error;
		}

	protected:
		MapVector m_timings;
		int m_cur_run_num;
		int m_total_runs;
		double m_l2_error;
		//debug counter, only outputted in debug mode
		std::map<int,int> m_count;


};

template < class CollectiveCommunication >
void Profiler::Output( CollectiveCommunication& comm, const int refineLevel, const long numDofs )
{
	const int numProce = comm.size();

	std::ostringstream filename;
	filename << "./data/p" << numProce << "_refinelvl_" << refineLevel << ".csv";
	filename.flush();

    if ( comm.rank() == 0 )
        std :: cout << "Profiling info in: " << ( filename.str() ).c_str() << std::endl;

#ifndef NDEBUG
    for (std::map<int,int>::const_iterator it = m_count.begin(); it != m_count.end(); ++it)
    {
        std::cout << "proc " << comm.rank() << " bId " << it->first << " count " << it->second << std::endl;
    }
#endif

    std::ofstream csv(( filename.str() ).c_str() );

    typedef std::map<std::string,long> AvgMap;
    AvgMap averages;
    for ( MapVector::const_iterator vit = m_timings.begin(); vit != m_timings.end(); ++vit )
	{
		for ( DataMap::const_iterator it = vit->begin(); it != vit->end(); ++it )
        {
            int diff = it->second.end - it->second.start;
            averages[it->first] +=diff;
        }
	}

//outputs column names
    csv << "refine;" << "processes; "<< "numDofs; " << "L1 error; ";
	for ( AvgMap::const_iterator it = averages.begin(); it != averages.end(); ++it )
	{
		csv << it->first << ";"  ;
	}
    csv << "Speedup (total); Speedup (ohne Solver)";

//outputs column values
	csv << std::endl
        << refineLevel << ";" << comm.size() << ";" << numDofs << ";"
        << comm.sum(m_l2_error)/double(m_total_runs*numProce) << ";";
	for ( AvgMap::const_iterator it = averages.begin(); it != averages.end(); ++it )
	{
		long clock_count = it->second;
		clock_count =  comm.sum( clock_count ) / double( CLOCKS_PER_SEC*0.001*numProce );
		csv << clock_count/double(m_total_runs) << ";" ;
	}
    csv << "=I$2/I2;" << "=SUM(E$2:G$2)/SUM(E2:G2)"  << std::endl;


	csv.close();

}

Profiler& profiler()
{
	static Profiler pf;
	return pf;

}

#endif // PROFILER_HH_INCLUDED
