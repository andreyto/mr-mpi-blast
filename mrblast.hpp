//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

////////////////////////////////////////////////////////////////////////////////
//
//  MR-MPI-BLAST: Parallelizing BLAST using MapReduce-MPI
//
//  Author: Seung-Jin Sul
//         (ssul@jcvi.org)
//
//  Last updated: 01/25/2012
//
////////////////////////////////////////////////////////////////////////////////

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

/// For importing search strategy
//#include <objects/blast/Blast4_request.hpp>
//#include <algo/blast/api/search_strategy.hpp>

/// For accessing fields in tabular output
//#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>

/// sequence::GetTitle() for get def line
#include <objmgr/util/sequence.hpp>

/// For getting align info and RecoverSearchStrategy
#include "blast_app_util.hpp"

/// For blastdbcmd, setting the dbsize
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>

/// For EScoreType
#include <objects/seqalign/Seq_align.hpp>

/// For sequence::GetLength()
#include <objmgr/util/seq_loc_util.hpp>

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

/// String trim
#include <boost/algorithm/string/trim.hpp>

/// ----------------------------------------------------------------------------
/// Settings from mrblast.ini file
/// ----------------------------------------------------------------------------
/// MR-MPI options
int g_verbosity;
int g_timer;                /// Log elapsed time for each mapreduce call
int g_memSize;              /// Page size (in Mbytes)
int g_outOfCore;

/// DB options
int g_nDbFiles;
string g_dbName;
string g_dbListFileName;
const string CONF_FILE_NAME = "mrblast.ini";
const int MAXPROCNAME = 80; /// For mpi proc name  
const int MAXSUBJID = 40;   /// For subject id string
const int MAXDBLENSTR = 40;   /// For database length string

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

typedef struct structQueryIndex {
    uint64_t qStart;                /// query start loc in orig query file
    uint32_t qLength;               /// quety len in bp
    uint64_t uniqQueryId;           /// unique query id (gi or serial num)
} structQueryIndex_t;

typedef struct structBlockBeginLoc {
    uint64_t blockBeginLoc;         /// query start loc of a work item block
    uint64_t qIdStart;              /// starting query id
} structBlockBeginLoc_t;

typedef struct structWorkItem {
    int      dbNo;                  /// db id
    uint64_t blockBegin;            /// query start loc of a work item block
    uint64_t blockEnd;              /// query end loc of a work item block
    uint64_t qIdStart;              /// starting query id
} structWorkItem_t;

vector<string> g_vecDbFile;
vector<structBlockBeginLoc_t> g_vecBlockBeginLoc; /// for generic version
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
int g_mapStyle;
int g_MPI_worldRank;        /// MPI rank
int g_MPI_nProcs;
char g_MPI_procName[MAXPROCNAME]; /// MPI procname
string g_outFilePrefix;     /// Prefix string for output file names
string g_indexFileName;
string g_queryFileName;

/// For nucl or prot DB setting
bool g_bIsProtein = false;

/// Add classifier fields in final BLAST output
bool g_bClassifier = false;

/// Output file to store BLAST hits
string g_hitFileName;
 
/// To emit MapReduce KV 
///
/// Fields
///    gi           : Query GI (GI of the origin sequence)
///    sId          : Subject Seq-id
///    identity     : BLAST percentage identity (Percentage of identical matches) (%)
///    alignLen     : alignment length
///    nMismatches  : number of mismatches
///    nGaps        : total number of gaps
///    queryStart   : Start of alignment in query
///    queryEnd     : End of alignment in query
///    subjectStart : Start of alignment in subject
///    subjectEnd   : End of alignment in subject
///    evalue       : Expect value
///    bitScore     : Bit score
///
///    identity     : Percentage identity for query (identity count/query length * 100, computed if classifier option is enabled)
///    coverage     : Percentage coverage for subject ((qend-qstart)/query length * 100, computed if classifier option is enabled)
///
typedef struct structBlResGeneric {
    char        subjectId[MAXSUBJID];
    double      identity;
    uint32_t    alignLen;
    uint32_t    nMismatches;
    uint32_t    nGaps;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    double      bitScore;
} structBlResGeneric_t;

/// +
/// percIdentity and percCoverage
typedef struct structBlResClassifier {
    char        subjectId[MAXSUBJID];
    double      identity;
    uint32_t    alignLen;
    uint32_t    nMismatches;
    uint32_t    nGaps;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    double      bitScore;
    double      percIdent;
    double      percCover; 
} structBlResClassifier_t;

/// To sort Blast hits by evalue & bitscore
typedef struct structEvalue {
    structBlResGeneric_t *pRec;
    double                eValue;
    double                bitScore;
    double                identity;
} structEValue_t;

/// +
/// percIdentity and percCoverage
typedef struct structEvalueClassifier {
    structBlResClassifier_t *pRec;
    double                   eValue;
    double                   bitScore;
    double                   identity;
} structEValueClassifier_t;

/// For saving hits in bin format
/// query id, subject id, % identity, alignment length, nMismatches, gap opens,
/// q. start, q. end, s. start, s. end, evalue, bit score
/// +
/// percIdentity and percCoverage
typedef struct structBlResToSaveHits {
    uint64_t    queryId;
    char        subjectId[MAXSUBJID];
    double      identity;
    uint32_t    alignLen;
    uint32_t    nMismatches;
    uint32_t    nGaps;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    double      bitScore;
} structBlResToSaveHits_t;

/// +
/// percIdentity and percCoverage
typedef struct structBlResToSaveHitsClassifier {
    uint64_t    queryId;
    char        subjectId[MAXSUBJID]; /// for the orig. defline of query
    double      identity;
    uint32_t    alignLen;
    uint32_t    nMismatches;
    uint32_t    nGaps;
    uint32_t    qStart;
    uint32_t    qEnd;
    uint32_t    sStart;
    uint32_t    sEnd;
    double      eValue;
    double      bitScore;
    double      percIdent;
    double      percCover;
} structBlResToSaveHitsClassifier_t;

/// For multiple iterations
int g_nIter = 1;   
int g_currIter = 0;
uint32_t nSubWorkItemSets = 0;
uint32_t nWorkItemsPerIter = 0;
uint32_t nRemains = 0;

/// For sorting hits by qid
uint32_t *g_vecNumHitsPerQid;
uint32_t *g_vecNumHitsPerQid2;

/// prototypes
void mr_mpi_blast(); 
void mr_map_run_blast(int itask, KeyValue *kv, void *ptr);
inline void mr_reduce_sort_and_save_generic(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline void mr_reduce_sort_and_save_classifier(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
inline int  mr_myhash(char* key, int len);
inline void mpi_collect_node_name(int rank, int nProcs, MPI_Comm mpiComm);
inline bool compare_evalue_generic(structEValue_t e1, structEValue_t e2);
inline bool compare_evalue_classifier(structEValueClassifier_t e1, structEValueClassifier_t e2);

#endif
