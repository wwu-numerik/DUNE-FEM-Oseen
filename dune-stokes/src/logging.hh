/**
 *  \file logging.hh
 *  \brief  logging
 **/
 #ifndef LOGGING_HH_INCLUDED
 #define LOGGING_HH_INCLUDED

 #include <fstream>
 #include <ostream>
 #include <sstream>
 #include <ctime>
 #include "stuff.hh"

class LogStream : virtual public std::ostream
{
    public:
        enum LogLevel{
          MIN = 1,
          ALL = 2,
          ERR = 4
        };

        LogStream( std::ofstream& logfile, LogLevel global_l, LogLevel local_l, bool log_to_file )
            :logfile_(logfile)
        {
            global_loglevel_ = global_l;
            this_loglevel_ = local_l;
            do_output = global_loglevel_ >= this_loglevel_;
            log_to_file_ = log_to_file;
        }

        ~LogStream()
        {
            //only flush, closing is done in parent
            logfile_.flush();
            std::cout.flush();
        }

        template < typename T >
        std::ostream& operator << ( T& input )
        {
            if ( do_output ) {
                if ( log_to_file_ ) {
                    logfile_ << input;
                }

                std::cout << input;
            }

            return std::cout;
        }


        std::ostream& operator<< (std::ostream& ( *pf )(std::ostream&))
        {
            if ( do_output ) {
                if ( log_to_file_ ) {
                    logfile_ << pf;
                }

                (this_loglevel_ == ERR ? std::cerr : std::cout ) << pf;
            }

            return (this_loglevel_ == ERR ? std::cerr : std::cout ) ;
        }

    private:
        bool log_to_file_,do_output;
        std::ofstream& logfile_;
        LogLevel global_loglevel_;
        LogLevel this_loglevel_;
};

class Logging
{
    public:


        Logging( )
            :   stream_min_(0),
                stream_max_(0),
                stream_err_(0)
        {
        }


        ~Logging()
        {
            Stuff::safe_delete( stream_max_ );
            Stuff::safe_delete( stream_min_ );
            Stuff::safe_delete( stream_err_ );

            if ( log_to_file_ ) {
                logfile_ << TimeString() << ": LOG END" << std::endl;
                logfile_.close();
            }
        }

        void Create (LogStream::LogLevel l, bool log_to_file = 0, std::string logfile = "log" )
        {
            loglevel_ = l;
            log_to_file_ = log_to_file;
            filename_ = logfile;
            if ( log_to_file_ ) {
                logfile_.open ( filename_.c_str() );
            }

            stream_err_ = new LogStream( logfile_ , l, LogStream::ERR , log_to_file );
            stream_min_ = new LogStream( logfile_ , l, LogStream::MIN , log_to_file );
            stream_max_ = new LogStream( logfile_ , l, LogStream::ALL , log_to_file );
        }

        LogStream& Min( )
        {
            return *stream_min_;
        }


        void Log( LogStream::LogLevel l, const std::string&  )
        {


        }

    private:
        bool log_to_file_;
        std::string filename_;
        std::ofstream logfile_;
        LogStream::LogLevel loglevel_;
        LogStream* stream_min_;
        LogStream* stream_max_;
        LogStream* stream_err_;

        std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }
};

static Logging& Logger ()
{
    static Logging log;
    return log;
}

#endif
