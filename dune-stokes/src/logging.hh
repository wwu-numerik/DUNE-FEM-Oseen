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
 #include <iomanip>
 #include "stuff.hh"
 #include "parametercontainer.hh"


/** \brief handles all logging
**/
class Logging
{
    public:

        enum LogFlags{
          LOG_ERR = 1,
          LOG_INFO = 2,
          LOG_DEBUG = 4,
          LOG_CONSOLE = 8,
          LOG_FILE = 16
        };

        class LogStream //: virtual public std::ostream
        {
            protected:
                LogFlags loglevel_;
                int& logflags_;
                std::stringstream buffer_;
                std::ofstream& logfile_;

            public:
                LogStream( LogFlags loglevel, int& logflags, std::ofstream& file )
                    : loglevel_(loglevel), logflags_(logflags),
                      logfile_(file) {}
                ~LogStream(){}
                template < typename T >
                LogStream& operator << ( T in )
                {
                    if ( logflags_ & loglevel_ )
                        buffer_ << in;
                    return *this;
                }

                //template < class Class >
                LogStream& operator << ( LogStream& ( *pf )(LogStream&) )
                {
                    if ( logflags_ & loglevel_ )
                        buffer_ << pf;
                    return *this;
                }

                LogStream& operator << ( std::ostream& ( *pf )(std::ostream &) )
                {
                    if ( logflags_ & loglevel_ ) {
                        if ( pf == (std::ostream& ( * )(std::ostream&))std::endl )
                        { //flush buffer into stream
                            if ( ( logflags_ & LOG_CONSOLE ) != 0 ) {
                                std::cout << buffer_.str() << std::endl;
                            }
                            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                                logfile_ << "\n" << TimeString()
                                         << buffer_.str() << std::endl;
                            }
                            buffer_.str("");// clear the buffer

                        }
                        else
                            buffer_ << pf;
                    }
                    return *this;
                }
        };

        Logging( )
        {
        }

        ~Logging()
        {
            Stuff::safe_delete( streammap_[LOG_INFO] );
            Stuff::safe_delete( streammap_[LOG_DEBUG] );
            Stuff::safe_delete( streammap_[LOG_ERR] );

            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_ << '\n' << TimeString() << ": LOG END" << std::endl;
                logfile_.close();
            }
        }


        /** \brief setup loglevel, logfilename
            \param logflags any OR'd combination of flags
            \param logfile filename for log, can contain paths, but creation will fail if dir is non-existant
        **/
        void Create (unsigned int logflags = LOG_FILE | LOG_CONSOLE | LOG_ERR, std::string logfile = "dune-stokes.log" )
        {
            logflags_ = logflags;
            filename_ = logfile;
            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_.open ( filename_.c_str() );
                assert( logfile_.is_open() );
            }
            flagmap_[LOG_ERR] = logflags;
            flagmap_[LOG_INFO] = logflags;
            flagmap_[LOG_DEBUG] = logflags;
            streammap_[LOG_ERR] = new LogStream( LOG_ERR, flagmap_[LOG_ERR], logfile_ );
            streammap_[LOG_DEBUG] = new LogStream( LOG_DEBUG, flagmap_[LOG_DEBUG], logfile_ );
            streammap_[LOG_INFO] = new LogStream( LOG_INFO, flagmap_[LOG_INFO], logfile_ );
        }

        void SetStreamFlags( LogFlags stream, int flags )
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            //this might result in logging to diff targtes, so we flush the current targets
            flagmap_[stream] = flags;
        }

         /** \name Log funcs for member-function pointers
         * \{
         */
        template < class Class >
        void LogDebug( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_DEBUG ) )
                Log( pf, c, LOG_DEBUG );
        }

        template < class Class >
        void LogInfo( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_INFO ) )
                Log( pf, c, LOG_INFO );
        }

        template < class Class >
        void LogErr( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_ERR ) )
                Log( pf, c, LOG_ERR  );
        }
        /** \}
        */

         /** \name Log funcs for basic types/classes
         * \{
         */
        template < class Class >
        void LogDebug( Class c )
        {
            if ( ( logflags_ & LOG_DEBUG ) )
                Log( c, LOG_DEBUG );
        }

        template < class Class >
        void LogInfo( Class c )
        {
            if ( ( logflags_ & LOG_INFO ) )
                Log( c, LOG_INFO );
        }

        template < class Class >
        void LogErr( Class c )
        {
            if ( ( logflags_ & LOG_ERR ) )
                Log( c, LOG_ERR );
        }

        LogStream& Err() { assert( streammap_[LOG_ERR] ); return *streammap_[LOG_ERR]; }
        LogStream& Info() { assert( streammap_[LOG_INFO] ); return *streammap_[LOG_INFO]; }
        LogStream& Dbg() { assert( streammap_[LOG_DEBUG] ); return *streammap_[LOG_DEBUG]; }
        /** \}
        */

        static std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }

    private:
        std::string filename_;
        std::ofstream logfile_;
        typedef std::map<LogFlags,int> FlagMap;
        FlagMap flagmap_;
        typedef std::map<LogFlags,LogStream*> StreamMap;
        StreamMap streammap_;
        int logflags_;

        template < class Class >
        void Log( void ( Class::*pf )(std::ostream&) , Class& c, LogFlags stream)
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            if ( ( flagmap_[stream] & LOG_CONSOLE ) != 0 )
                (c.*pf)( std::cout );
            if ( ( flagmap_[stream] & LOG_FILE ) != 0 )
                (c.*pf)( logfile_ );
        }

        template < class Class >
        void Log( Class c, LogFlags stream )
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            if ( ( flagmap_[stream] & LOG_CONSOLE ) != 0 )
                std::cout << c;
            if ( ( flagmap_[stream] & LOG_FILE ) != 0 )
                logfile_ << c;
        }
};

///global Logging instance
Logging& Logger ()
{
    static Logging log;
    return log;
}

#endif
