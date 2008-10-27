/**
 *  \file logging.hh
 *  \brief  logging
 **/
 #include <fstream>
 #include <ostream>
 #include <sstream>
 #include <ctime>

class LogStream : public std::ostream
{
    public:
        enum LogLevel{
          MIN,
          ALL,
          ERR
        };

        LogStream( std::ofstream& logfile, LogLevel l, bool log_to_file )
            :logfile_(logfile)
        {
            loglevel_ = l;
            log_to_file_ = log_to_file;
        }

        template < typename T >
        std::ostream& operator << ( T& input )
        {
            return buffer_ << input;
        }


        std::ostream& operator << ( void*& input )
        {
            if ( log_to_file_ )
            {
                logfile_ << buffer_ << std::endl;
            }

            buffer_.clear();
            return buffer_;
        }

    private:
        bool log_to_file_;
        std::ofstream& logfile_;
        LogLevel loglevel_;
        std::stringstream buffer_;

};

class Logging
{
    public:


        Logging( )
            :   stream_min_(0),
                stream_max_(0),
                stream_err_(0)
         { }
        ~Logging()
        {
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
        }

        void Min( const std::string& msg )
        {
            Log( MIN, msg );
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



// function.print( Logging::Min )
