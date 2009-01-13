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
          LOG_NONE = 1,
          LOG_ERR = 2,
          LOG_INFO = 4,
          LOG_DEBUG = 8,
          LOG_CONSOLE = 16,
          LOG_FILE = 32
        };

        class LogStream //: virtual public std::ostream
        {
            protected:
                LogFlags loglevel_;
                int& logflags_;
                std::stringstream buffer_;
                std::ofstream& logfile_;
                std::ofstream& logfileWoTime_;

            public:
                LogStream( LogFlags loglevel, int& logflags, std::ofstream& file, std::ofstream& fileWoTime )
                    : loglevel_(loglevel), logflags_(logflags),
                      logfile_(file), logfileWoTime_( fileWoTime ) {}
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
                                logfileWoTime_ << buffer_.str() << std::endl;
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
            streamIDs_.push_back( LOG_ERR );
            streamIDs_.push_back( LOG_DEBUG );
            streamIDs_.push_back( LOG_INFO );
        }

        ~Logging()
        {
            IdVecCIter it = streamIDs_.end();
            for ( ; it != streamIDs_.begin(); --it ) {
                Stuff::safe_delete( streammap_[*it] );
            }

            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_ << '\n' << TimeString() << ": LOG END" << std::endl;
                logfileWoTime_ << std::endl;
                logfile_.close();
                logfileWoTime_.close();
            }
        }


        /** \brief setup loglevel, logfilename
            \param logflags any OR'd combination of flags
            \param logfile filename for log, can contain paths, but creation will fail if dir is non-existant
        **/
        void Create (unsigned int logflags = LOG_FILE | LOG_CONSOLE | LOG_ERR, std::string logfile = "dune_stokes" )
        {
            logflags_ = logflags;
            filename_ = logfile + "_time.log";
            filenameWoTime_ = logfile + ".log";
            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_.open ( filename_.c_str() );
                assert( logfile_.is_open() );
                logfileWoTime_.open ( filenameWoTime_.c_str() );
                assert( logfileWoTime_.is_open() );
            }
            IdVecCIter it = streamIDs_.begin();
            for ( ; it != streamIDs_.end(); ++it ) {
                flagmap_[*it] = logflags;
                streammap_[*it] = new LogStream( *it, flagmap_[*it], logfile_, logfileWoTime_ );
            }
        }

        void SetStreamFlags( LogFlags stream, int flags )
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            //this might result in logging to diff targtes, so we flush the current targets
            flagmap_[stream] = flags;
        }

        int GetStreamFlags( LogFlags stream )
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            return flagmap_[stream];
        }

         /** \name Log funcs for member-function pointers
         * \{
         */
        template < typename Pointer, class Class >
        void LogDebug( Pointer pf, Class& c )
        {
            if ( ( logflags_ & LOG_DEBUG ) )
                Log( pf, c, LOG_DEBUG );
        }

        template < class Class,typename Pointer >
        void LogInfo( Pointer pf , Class& c )
        {
            if ( ( logflags_ & LOG_INFO ) )
                Log( pf, c, LOG_INFO );
        }

        template < typename Pointer, class Class >
        void LogErr( Pointer pf, Class& c )
        {
            if ( ( logflags_ & LOG_ERR ) )
                Log( pf, c, LOG_ERR  );
        }

        template < typename Pointer, class Class >
        void Log( void ( Class::*pf )(std::ostream&) const, Class& c, LogFlags stream)
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            if ( ( flagmap_[stream] & LOG_CONSOLE ) != 0 )
                (c.*pf)( std::cout );
            if ( ( flagmap_[stream] & LOG_FILE ) != 0 )
                (c.*pf)( logfile_ );
        }

        template < class Class,typename Pointer >
        void Log( Pointer pf, Class& c, LogFlags stream)
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            if ( ( flagmap_[stream] & LOG_CONSOLE ) != 0 )
                (c.*pf)( std::cout );
            if ( ( flagmap_[stream] & LOG_FILE ) != 0 ) {
                (c.*pf)( logfile_ );
                (c.*pf)( logfileWoTime_ );
            }
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

        template < class Class >
        void Log( Class c, LogFlags stream )
        {
            assert( stream & ( LOG_ERR | LOG_INFO | LOG_DEBUG ) );
            if ( ( flagmap_[stream] & LOG_CONSOLE ) != 0 )
                std::cout << c;
            if ( ( flagmap_[stream] & LOG_FILE ) != 0 ) {
                logfile_ << c;
                logfileWoTime_ << c;
            }
        }
        /** \}
        */

        LogStream& GetStream(int stream) { assert( streammap_[(LogFlags)stream] ); return *streammap_[(LogFlags)stream]; }
        LogStream& Err() { return GetStream( LOG_ERR ); }
        LogStream& Info() { return GetStream( LOG_INFO ); }
        LogStream& Dbg() { return GetStream( LOG_DEBUG ); }

        static std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }

        int AddStream( int flags )
        {
//            assert( streamIDs_.find( streamID ) == streamIDs_.end() );
            static int streamID_int = 16;
            streamID_int <<= 2;
            LogFlags streamID = (LogFlags) streamID_int;
            streamIDs_.push_back( streamID );
            flagmap_[streamID] = flags | streamID;
            streammap_[streamID] = new LogStream( streamID, flagmap_[streamID], logfile_, logfileWoTime_ );
            return streamID_int;
        }

    private:
        std::string filename_;
        std::string filenameWoTime_;
        std::ofstream logfile_;
        std::ofstream logfileWoTime_;
        typedef std::map<LogFlags,int> FlagMap;
        FlagMap flagmap_;
        typedef std::map<LogFlags,LogStream*> StreamMap;
        StreamMap streammap_;
        typedef std::vector<LogFlags> IdVec;
        typedef std::vector<LogFlags>::const_iterator IdVecCIter;
        IdVec streamIDs_;
        int logflags_;
};

///global Logging instance
Logging& Logger ()
{
    static Logging log;
    return log;
}

#endif
