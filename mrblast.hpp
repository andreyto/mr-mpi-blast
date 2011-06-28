#ifndef MRBLAST_H
#define MRBLAST_H

/// For typedef unsigned long long int uint32_t
#include <stdint.h>

/// MPI AND MR-MPI
#include "mpi.h"
#include "mrmpi/mapreduce.h"
#include "mrmpi/keyvalue.h"

/// For syncronized timing
#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL 1
#endif

using namespace MAPREDUCE_NS;
using namespace std;

/// ----------------------------------------------------------------------------
/// NCBI C++ Toolkit
/// ----------------------------------------------------------------------------
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

/// Import search strategy
#include <objects/blast/Blast4_request.hpp>
#include <algo/blast/api/search_strategy.hpp>

/// Tabular
#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>

/// sequence::GetTitle() for get def line
#include <objmgr/util/sequence.hpp>

/// For align info
#include "blast_app_util.hpp"

USING_NCBI_SCOPE;
USING_SCOPE (blast);
USING_SCOPE (align_format);

const int QUERY = 0;        /// To retireve query info from CSeq_align
const int SUBJECT = 1;      /// To retireve suject info from CSeq_align

/// For NcbiApp
CRef<CBlastOptionsHandle> g_opts_hndl;
CRef<CBlastAppArgs> g_cmdLineArgs;

/// Blast target DB setting
CRef<CSearchDatabase> g_searchDatabase;
string g_prevDbName = "";
    
/// ----------------------------------------------------------------------------
/// Boost lib
/// ----------------------------------------------------------------------------
/// For tokenizer
#include <boost/algorithm/string.hpp>

/// For str -> int, int -> str
#include <boost/lexical_cast.hpp>

/// For processing command line arguments
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// For processing configuration file
#include <boost/program_options/detail/config_file.hpp>
namespace pod = boost::program_options::detail;

/// For Boost memory mapped file
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/operations.hpp> /// for real file size
boost::iostreams::mapped_file_source g_memmapQueryFile; /// Read-only Boost mmap query file
uintmax_t g_realFileSize = 0;

/// ----------------------------------------------------------------------------
/// Settings from mrblast.ini conf file
/// ----------------------------------------------------------------------------
/// MR-MPI options
int g_verbosity;
int g_timer;                /// log elapsed time for each mapreduce call
int g_memSize;              /// page size (in Mbytes)
int g_outOfCore;

/// DB options
string g_dbFileName;
string g_configFileName = "mrblast.ini";
int g_numDbFiles;
const int MAXSTR = 80;      /// For mpi proc name 

/// Filtering
float g_IDENT_CUTOFF = 0.0; /// Doug's identity for filtering
float g_COVER_CUTOFF = 0.0; /// Doug's coverage for filtering

/// ----------------------------------------------------------------------------
/// Log
/// ----------------------------------------------------------------------------
#include <sys/time.h>
#include <sys/resource.h>
//#define NDEBUG 1
int g_bLogEnabled = 0;
int g_timingEnabled = 0;
int g_optDumpEnabled = 0;   /// For dumping Blast opts_hndl out
int g_mapCallNo = 0;        /// Serial no of map calls in each rank
string g_logFileName;
string g_logMsg;
ofstream g_logFileStream;
static ostream *g_pLogStream = NULL;
#define LOG (*g_pLogStream)

/// ----------------------------------------------------------------------------
/// New scheduler
/// ----------------------------------------------------------------------------
typedef multimap<string, int> multimapSI_t;
typedef multimap<string, vector<int> > multimapSVI_t;
typedef pair<string, vector<int> > pairSVI_t;
typedef multimap<int, int> multimapII_t;
typedef map<int, string> mapIS_t;
typedef map<string, string> mapSS_t;

typedef struct structIndex {
    uintmax_t qStart;
    uint32_t qLength;
} structQueryIndex_t;

typedef struct structWorkItem {
    uint32_t dbNo;
    uintmax_t blockBegin;
    uintmax_t blockEnd;
} structWorkItem_t;

vector<string> g_vecDbFile;
vector<uintmax_t> g_vecBlockBeginLoc;
vector<structQueryIndex_t> g_vecQueryIndex;
vector<structWorkItem_t> g_vecWorkItem;
uint32_t g_numQueries;
uint32_t g_numQueryBlocks;
uint32_t g_numWorkItems;
uint32_t g_blockSize;                         /// block size in base-pair
multimapSI_t g_multimapProcNameRank; /// dict of rank by proc name
mapIS_t g_mapRankProcName;           /// dict of proc name by rank 
 
/// ----------------------------------------------------------------------------
/// Misc.
/// ----------------------------------------------------------------------------
string g_outFilePrefix;     /// Prefix string for output file names
string g_indexFileName;
string g_queryFileName;
int g_mapStyle;
int g_MPI_worldRank;         /// MPI rank
char g_MPI_procName[MAXSTR]; /// MPI procname
int g_MPI_numProcs;

/// For nucl or prot DB setting
bool g_bIsProtein = false;

/// Output file to store BLAST hits
string g_hitFileName = "";
 
/// To pass Blast hits following outfmt=6 format.
/// subject id, % identity, alignment length, nMismatches, gap opens,
/// q. start, q. end, s. start, s. end, evalue, bit score
typedef struct structBlRes {
    uint32_t    subjectId;    
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       doug_identity;
    float       doug_coverage;
} structBlRes_t;

/// For saving hits in bin format
typedef struct structBlRes2 {
    uint32_t    qId;
    uint32_t    subjectId;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       doug_identity;
    float       doug_coverage;
} structBlRes2_t;

/// To sort Blast hits by evalue & bitscore
typedef struct structEvalue {
    structBlRes_t *pRec;
    uint32_t       subjectId;
    double         eValue;
    float          bitScore;
} structEValue_t;

/// For multiple iterations
uint16_t g_numIter = 1;   
uint16_t g_currIter = 0;
uint32_t nSubWorkItemSets = 0;
uint32_t nWorkItemsPerIter = 0;
uint32_t nRemains = 0;


#endif
