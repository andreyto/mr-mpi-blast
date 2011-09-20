//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#ifndef MRBLAST_HPP
#define MRBLAST_HPP

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
CRef<CBlastOptionsHandle> g_optsHndl;
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
uint64_t g_realFileSize = 0;

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
string CONF_FILE_NAME = "mrblast.ini";
int g_nDbFiles;
const int MAXSTR = 80;      /// For mpi proc name 

/// Filtering
float g_identCutoff = 0.0; /// Doug's identity for filtering
float g_coverCutoff = 0.0; /// Doug's coverage for filtering

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
    uint64_t qStart;
    uint32_t qLength;
    uint64_t uniqQueryId;
} structQueryIndex_t;

typedef struct structWorkItem {
    int dbNo;
    uint64_t blockBegin;
    uint64_t blockEnd;
} structWorkItem_t;

vector<string> g_vecDbFile;
vector<uint64_t> g_vecBlockBeginLoc;
vector<structQueryIndex_t> g_vecQueryIndex;
vector<structWorkItem_t> g_vecWorkItem;
uint32_t g_nQueries;
uint32_t g_nQueryBlocks;
uint32_t g_nWorkItems;
uint32_t g_blockSize;                /// block size in base-pair
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
int g_MPI_nProcs;

/// For nucl or prot DB setting
bool g_bIsProtein = false;

/// Output file to store BLAST hits
string g_hitFileName;
 
/// To pass Blast hits following outfmt=6 format.
/// subject id, % identity, alignment length, nMismatches, gap opens,
/// q. start, q. end, s. start, s. end, evalue, bit score
typedef struct structBlRes {
    uint64_t    subjectId;    
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       identity;
    float       coverage;
} structBlRes_t;

typedef struct structBlResMason {
    uint64_t    subjectId;    
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    float       identity;
    float       coverage;
} structBlResMason_t;

/// For saving hits in bin format
typedef struct structBlResToSaveHits {
    uint64_t    gi;
    uint64_t    subjectId;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       identity;
    float       coverage;
} structBlResToSaveHits_t;

/// KV key = gi_readid_cutX_cutY_F/R_P0/P1
typedef struct structBlResToSaveHitsMason {
    uint64_t    gi;
    uint64_t    readId;
    char        readStrand;
    char        readPairId;    
    uint64_t    subjectId;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    float       bitScore;
    uint32_t    upperStart;
    uint32_t    upperEnd;
    float       identity;
    float       coverage;    
} structBlResToSaveHitsMason_t;

/// To sort Blast hits by evalue & bitscore
typedef struct structEvalue {
    structBlRes_t *pRec;
    uint64_t       subjectId;
    double         eValue;
    float          bitScore;
} structEValue_t;

/// For multiple iterations
int g_nIter = 1;   
int g_currIter = 0;
uint32_t nSubWorkItemSets = 0;
uint32_t nWorkItemsPerIter = 0;
uint32_t nRemains = 0;


void mrmpi_blast(); 
void mr_map_run_blast(int itask, KeyValue *kv, void *ptr);
inline void mpi_collect_node_name(int rank, int nProcs, MPI_Comm mpiComm);
inline void mr_reduce_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline void mr_reduce_save_hits(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline bool compare_evalue(structEValue_t e1, structEValue_t e2);


#endif
