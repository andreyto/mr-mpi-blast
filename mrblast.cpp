////////////////////////////////////////////////////////////////////////
//
//  MPI-MR-BLAST: Parallelizing BLAST on MR-MPI
//
//  Author: Seung-Jin Sul
//          (ssul@jcvi.org)
//
////////////////////////////////////////////////////////////////////////

/*
 *     This program is free software; you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation; either version 2 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program; if not, write to the Free Software
 *     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA 02110-1301, USA.
 */

/*
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 */

/// MPI AND MR-MPI
#include "mpi.h"
#include "mrmpi/mapreduce.h"
#include "mrmpi/keyvalue.h"

using namespace MAPREDUCE_NS;
using namespace std;

/// NCBI
#include <ncbi_pch.hpp>

#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>
#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

/// ASN & split
#include <boost/algorithm/string.hpp>

/// For args processing for Blast
#include <algo/blast/blastinput/blastn_args.hpp>

/// For typedef unsigned long long int uint64_t
#include <stdint.h>

/// For queue
#include <queue>

/// str2int, int2str
#include <boost/lexical_cast.hpp>

/// Processing command line arguments
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// Configuration file processing
#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
namespace pod = boost::program_options::detail;

/// Import search strategy
#include <objects/blast/Blast4_request.hpp>
#include <algo/blast/api/search_strategy.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);

/// Etc
const int MAXSTR = 256;
const int QUERY = 0;
const int SUBJECT = 1;

/// ----------------------------------------------------------------------------
/// Set from mrblast.ini conf file
/// ----------------------------------------------------------------------------
unsigned int NMAXQUERY; /// num query to accumulate from qeury file
int EXCLUSIONTHRESHOLD; /// Exclusion threshold = 100bp

/// Blast options
double EVALUE;
string BLASTOPTS; /// To pass Blast user options 

/// DB options
string DBCHUNKLISTFILE;
string EXCLUSIONHISTFILE;
int NTOTALDBCHUNKS;
int NITERATION;

/// Log
int LOGORNOT = 0;
int LOGTOFILE = 0;
ostream* LOGSTREAM = NULL;
#define LOG (*LOGSTREAM)
string LOGMSG;
string LOGFILENAME;
/// ----------------------------------------------------------------------------

/// Blast target DB setting
static CSearchDatabase *pTargetDb = 0;
string prevDbChunkName;

/// Import search strategy
string STRATEGYFILENAME;

/// Misc.
unsigned int MYID;
char* PNAME;
bool EXCLUSIONORNOT;
string OUTPREFIX;
string INDEXFILENAME;
string QUERYFILENAME;
bool OPTDUMPED = false; /// For loggin Blast opts
#define NDEBUG 1

/// For timing
double prog_start;
#define MPI_WTIME_IS_GLOBAL 1

/// To pass Blast hits, outfmt=6
/// subject id, % identity, alignment length, mismatches, gap opens, 
/// q. start, q. end, s. start, s. end, evalue, bit score
typedef struct blastres {
    unsigned long long subjectid;
    double identity;
    int alignlen; 
    int mismatches;
    int gapopens;
    int qstart; 
    int qend;
    int sstart;
    int send;
    double evalue; 
    int bitscore;
} BLASTRES;
int BLASTRESSZ = sizeof(BLASTRES);

/// To sort Blast hits by evalue
typedef struct structEvalue {
    BLASTRES* p;
    double evalue;
} STRUCTEVALUE;
string WORKEROUTPUTFILENAME = "";

/// To pass NTOTALDBCHUNKS and NCOREPERNODE to map() for custom scheduler
/// NMAXTRIAL = max num of trial to find node number which 
///             has the DB chunk of the work item.
typedef struct gf {
    int NTOTALDBCHUNKS;
    int NCOREPERNODE;
    int NMAXTRIAL;
} GF;
GF GF1;
 
/// MR-MPI CALLS
void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
void    mr_run_blast(int itask, char *file, KeyValue *kv, void *ptr);                   
void    mr_sort_multivalues_by_evalue(char *key, int keybytes, char *multivalue, 
            int nvalues, int *valuebytes, KeyValue *kv, void *ptr);                         

/// Check hit exclusion
bool    check_exclusion(string qGi, string sGi, uint64_t qCutLocStart, 
            uint64_t qCutLocEnd, uint64_t sStart, uint64_t sEnd, int threshold);

/* -------------------------------------------------------------------------- */
int main(int argc, char **argv)
/* -------------------------------------------------------------------------- */
{
    ///    
    /// Read conf file, mrblast.ini and set parameters
    ///
    ifstream config("mrblast.ini", ios::in);
    if (!config) {
        cerr << "ERROR: configuration file, mrblast.ini, not found" << endl;
        return 1;
    }
    
    std::set<std::string> options;
    std::map<std::string, std::string> parameters;
    options.insert("*");
    
    try {      
        for (pod::config_file_iterator i(config, options), e ; i != e; ++i) {
            parameters[i->string_key] = i->value[0];
        }
        DBCHUNKLISTFILE = parameters["DBCHUNKLISTFILE"];
        EXCLUSIONHISTFILE = parameters["EXCLUSIONHISTFILE"];
        LOGFILENAME = parameters["LOGFILENAME"];
        BLASTOPTS = parameters["BLASTOPTS"];
        try {   
            NMAXQUERY = boost::lexical_cast<unsigned int>(
                parameters["NMAXQUERY"]);
            NTOTALDBCHUNKS = boost::lexical_cast<unsigned int>(
                parameters["NTOTALDBCHUNKS"]);
            GF1.NTOTALDBCHUNKS = NTOTALDBCHUNKS;
            GF1.NCOREPERNODE = boost::lexical_cast<unsigned int>(
                parameters["NCOREPERNODE"]);
            GF1.NMAXTRIAL= boost::lexical_cast<int>(
                parameters["NMAXTRIAL"]);           
            EVALUE = boost::lexical_cast<double>(parameters["EVALUE"]);
            EXCLUSIONTHRESHOLD = boost::lexical_cast<int>(
                parameters["EXCLUSIONTHRESHOLD"]);     
            LOGORNOT = boost::lexical_cast<int>(parameters["LOGORNOT"]);
            LOGTOFILE = boost::lexical_cast<int>(parameters["LOGTOFILE"]);
        }
        catch(const boost::bad_lexical_cast &) {
            cerr << "Exception: bad_lexical_cast" << endl;
        }
    }
    catch(exception& e) {
        cerr<< "Exception: " << e.what() << endl;
    }
        
    /// 
    /// nquery: number of queries for a work item. If there are q queries 
    /// in the input file, (q/nquery) blocks of queries are generated. In turn, 
    /// If there are p DB partitions, there are p*(q/nquery) work items.
    /// nquery=0 means there is only one block of queries and therefore p work
    /// items will be generated.
    ///
    /// NOTE: query bp size should be considered for nquery later.
    ///
    /// iteration: number of iterations. p*(q/nquery) work items are divided 
    /// into p*(q/nquery) / iteration files and each will be processed one by 
    /// one.
    ///
    /// If mapstyle = 2, to use the default master/slave mode in MR-MPI.
    /// If mapstyle = 3, to use the modified master/slave mode for controlling 
    /// DB chunk locality.
    ///
    int MAPSTYLE = 0;
    unsigned int NQUERYPERWORKITEM;
    
    po::options_description generalDesc("General options");
    generalDesc.add_options() 
        ("help", "print help message")        
        ("query-file,i", po::value<string>(), "set input query file")
        ("index-file,d", po::value<string>(), "set input index file")
        ("import_search_strategy,s", po::value<string>(), 
            "set search strategy file")
    ;
    
    po::options_description OptionalDesc("Optional");
    OptionalDesc.add_options() 
        ("block-size,b", po::value<unsigned int>(&NQUERYPERWORKITEM)->default_value(0), 
            "set the number/size of queries per work item (default=0, which means all queries)")
        ("iteration,n", po::value<int>(&NITERATION)->default_value(1), 
            "set the number of iterations")            
        ("output-prefix,o", po::value<string>(&OUTPREFIX)->default_value("output"), 
            "set output prefix for output file names (default=output)")
        ("map-style,m", po::value<int>(&MAPSTYLE)->default_value(2), 
            "set MR-MPI mapstyle: 2=master/slave, 3=new scheduler")
        ("self-exclusion,x", po::value<bool>(&EXCLUSIONORNOT)->default_value(false), 
            "enable self exclusion")
    ;
    
    po::options_description allDesc("Allowed options");
    allDesc.add(generalDesc).add(OptionalDesc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, allDesc), vm);
    po::notify(vm); 
        
    if (argc < 2 || (!strcmp(argv[1], "-?") || !strcmp(argv[1], "--?") 
        || !strcmp(argv[1], "/?") || !strcmp(argv[1], "/h") 
        || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--h") 
        || !strcmp(argv[1], "--help") || !strcmp(argv[1], "/help") 
        || !strcmp(argv[1], "-help")  || !strcmp(argv[1], "help") )) {
        cout << "MR-MPI Blast\n"
             << "Author: Seung-Jin Sul (ssul@jcvi.org)\n\n"
             << allDesc;
        return 1;
    }
    else {
        if (vm.count("query-file")) {
            QUERYFILENAME = vm["query-file"].as<string>();
        }
        else {
            cerr << "ERROR: query file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("index-file")) {
            INDEXFILENAME = vm["index-file"].as<string>();
        }
        else {
            cerr << "ERROR: index file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        if (vm.count("import_search_strategy")) {
            STRATEGYFILENAME = vm["import_search_strategy"].as<string>();
        }
        else {
            cerr << "ERROR: index file was not set.\n\n"; 
            cout << allDesc;
            return 1;
        }
        
        /// 
        /// Read OPTIONAL arguments
        ///
        if (vm.count("output-prefix")) 
            OUTPREFIX = vm["output-prefix"].as<string>();
        if (vm.count("block-size")) 
            NQUERYPERWORKITEM = vm["block-size"].as<unsigned int>();
        if (vm.count("iteration")) 
            NITERATION = vm["iteration"].as<int>();
        if (vm.count("map-style")) 
            MAPSTYLE = vm["map-style"].as<int>();
        if (vm.count("self-exclusion")) 
            EXCLUSIONORNOT = vm["self-exclusion"].as<bool>();
    }
        
    ///
    /// Log stream setting
    ///
    ofstream logFile;
    if (LOGORNOT) {
        if (LOGTOFILE) {
            LOGFILENAME = OUTPREFIX + LOGFILENAME;
            logFile.open(LOGFILENAME.c_str(), ios::out);
            LOGSTREAM = &logFile;
        } 
        else LOGSTREAM = &cout;        
    }  
    
    ///
    /// MPI setup
    ///
    char MPI_procName[MAXSTR];
    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    
    LOGMSG = "[LOG] Rank:" + boost::lexical_cast<string>(MPI_myId) + " ";
    LOG << LOGMSG << "proc name = " << MPI_procName << endl;
    
    MPI_Barrier(MPI_COMM_WORLD);  
    double profile_time = MPI_Wtime();
    
    /// 
    /// MR-MPI init
    ///
    LOG << LOGMSG << "MR-MPI Init starts.\n";
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    
    /*
    * mapstyle = 0 (chunk) or 1 (stride) or 2 (master/slave)
    * all2all = 0 (irregular communication) or 1 (use MPI_Alltoallv)
    * verbosity = 0 (none) or 1 (summary) or 2 (histogrammed)
    * timer = 0 (none) or 1 (summary) or 2 (histogrammed)
    * memsize = N = number of Mbytes per page of memory
    * minpage = N = # of pages to pre-allocate per processor
    * maxpage = N = max # of pages allocatable per processor, 0 = no limit
    * keyalign = N = byte-alignment of keys
    * valuealign = N = byte-alignment of values
    * fpath = string 
    */
    mr2->verbosity = 0;
    mr2->timer = 0;
    mr2->mapstyle = MAPSTYLE;  /// master/slave mode=2, custom scheduler=3
    MPI_Barrier(MPI_COMM_WORLD);
    prog_start = MPI_Wtime();
    MYID = MPI_myId;
    PNAME = MPI_procName;
        
    ///
    /// Make a file for map() which contains a list of file
    /// names. A file name is a form of "queryFile,dbChunkName"
    /// In map(), the file name is splitted into "queryFile" and
    /// "dbChunkName". We've got 109 DB chunks, FYI.
    ///
    vector<string> vWorkItems;
    vector<string> vDbChunkNames;
    unsigned int nDbChunks;
    unsigned int nWorkItems;
    uint64_t nQueryBlocks = 0;
    
    double master_init_time = MPI_Wtime();
    LOG << LOGMSG << "Master's init work starts." << endl;
           
    if (MPI_myId == 0) {
        ///
        /// Read DB partition name list file
        ///
        string line;
        ifstream dbChunkNameFile(DBCHUNKLISTFILE.c_str(), ios::in);
        if (!dbChunkNameFile.is_open()) {
            cerr << "ERROR: dbchunks.txt open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        while (!getline(dbChunkNameFile, line).eof()) {
            vDbChunkNames.push_back(line);
        }
        dbChunkNameFile.close();        
        nDbChunks = vDbChunkNames.size();
        
        ///
        /// Read seqeunce line index file 
        /// Make work items from vDbChunkNames and vQueryFileNames
        /// a work item = <(begin_offset, end_offset), DBChunkName>
        /// *_offset mean the line index of the sequences in the 
        /// original query sequence file.
        ///
        ifstream indexFile(INDEXFILENAME.c_str(), ios::in);
        if (!indexFile.is_open()) {
            cerr << "ERROR: indexFile open error.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        string beginOffset;
        string endOffset;
        if (NQUERYPERWORKITEM) {
            unsigned int n = 0;
            while (!getline(indexFile, line).eof()) {
                vector<string> vOffsets;
                boost::split(vOffsets, line, boost::is_any_of(","));
                            
                if (n == 0) beginOffset = vOffsets[0];
                endOffset = vOffsets[1];
                
                if (n == NQUERYPERWORKITEM - 1) {
                    endOffset = vOffsets[1];
                    for (uint64_t j = 0; j < nDbChunks; ++j) { 
                        vWorkItems.push_back(beginOffset + "," + endOffset + ","
                            + vDbChunkNames[j]);
                    }
                    nQueryBlocks++;
                    n = 0;
                }      
                else n++;
            }        
            if (n != 0) {
                for (uint64_t j = 0; j < nDbChunks; ++j) { 
                    vWorkItems.push_back(beginOffset + "," + endOffset + ","
                        + vDbChunkNames[j]);
                }
                nQueryBlocks++;
            }
        }
        else {
            /// Use all
            getline(indexFile, line);
            vector<string> vOffsets;
            boost::split(vOffsets, line, boost::is_any_of(","));
            beginOffset = vOffsets[0];
            string temp;
            do { temp = line; } while (!getline(indexFile, line).eof());
            boost::split(vOffsets, temp, boost::is_any_of(","));   
            endOffset = vOffsets[1];   
            for (uint64_t j = 0; j < nDbChunks; ++j) { 
                vWorkItems.push_back(beginOffset + "," + endOffset + ","
                    + vDbChunkNames[j]);
            }       
            nQueryBlocks = 1;      
        }        
        indexFile.close();
        nWorkItems = vWorkItems.size();
        cout << "INFO: number of query blocks = " << nQueryBlocks << endl;     
        cout << "INFO: number of DB partitions = " << nDbChunks << endl;     
        cout << "INFO: number of work items = " << nWorkItems << endl;
 
        vDbChunkNames.clear();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);  
    
    ///
    /// Split work item master file into sub master files and process the 
    /// work items iteratively while saving the results at each iteration.
    /// eq) if the total number of query files x and NITERATION = x/2, 
    /// the master file is divided into 2 sub master files and processed 
    /// individually.
    ///
    uint64_t nSubWorkItemFiles = 0;
    uint64_t nWorkItemsPerFile = 0;
    if (MYID == 0) { 

        
        uint64_t nRemains = 0;
        if (NITERATION > 1) {
            nSubWorkItemFiles = NITERATION;
            nWorkItemsPerFile = (nQueryBlocks / NITERATION) * nDbChunks;            
            nRemains = (nQueryBlocks % NITERATION) * nDbChunks;
        }
        else { /// NITERATION = 1
            nWorkItemsPerFile = nWorkItems;
            nRemains = 0;
            nSubWorkItemFiles = 1;
        }
                
        uint64_t k = 0;
        uint64_t i = 0;
        
        for (; i < nSubWorkItemFiles; ++i) {
            string newFileName = OUTPREFIX + "-workitems-" 
                + boost::lexical_cast<string>(i) + ".txt";
            ofstream workItemFile(newFileName.c_str());
            cout << "INFO: Work item file name (" << i << ") = " 
                 << newFileName << endl;
            if (workItemFile.is_open()) {
                for (uint64_t j = 0; j < nWorkItemsPerFile; ++j, ++k) {
                    workItemFile << vWorkItems[k] << endl;
                }
                workItemFile.close();
            }
            else {
                cerr << "ERROR: work item file open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        if (nRemains) {
            string newFileName = OUTPREFIX + "-workitems-" 
                + boost::lexical_cast<string>(i) + ".txt";
            ofstream workItemFile(newFileName.c_str());
            cout << "INFO: Work item file name (remains) = " 
                 << newFileName << endl;
            if (workItemFile.is_open()) {
                for (uint64_t j = 0; j < nRemains; ++j, ++k) {
                    workItemFile << vWorkItems[k] << endl;
                }
                workItemFile.close();
                nSubWorkItemFiles++; /// for one remains file
            }
            else {
                cerr << "ERROR: work item file open error.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        cout << "INFO: Total number of sub work item files = " 
             << nSubWorkItemFiles << endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD); 
    LOG << LOGMSG << "Master's init work ends." 
        << "\t" << MPI_Wtime() - master_init_time << endl;
    
    ///
    /// Iteratively call blast and save results for nSubWorkItemFiles
    /// times.
    ///

    /// Broadcast the number of sub work item files
    MPI_Bcast(&nSubWorkItemFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    uint64_t nvecRes;
        
    for (uint64_t n = 0; n < nSubWorkItemFiles; ++n) {
        vector<string> vSubWorkitems;
        string subMasterFileName 
            = OUTPREFIX + "-workitems-" + boost::lexical_cast<string>(n) 
            + ".txt";
   
        ///
        /// map, collate reduce
        ///
        double map_time = MPI_Wtime();
        LOG << LOGMSG << "map() starts." << endl;
        nvecRes = mr2->map((char*)subMasterFileName.c_str(), &mr_run_blast, &GF1);
        LOG << LOGMSG << "map() ends." << "\t" <<  MPI_Wtime() - map_time << endl;
                            
        WORKEROUTPUTFILENAME 
            = OUTPREFIX + "-" + boost::lexical_cast<string>(n) + "-" 
            + boost::lexical_cast<string>(MPI_myId) + ".txt";
        
        double collate_time = MPI_Wtime();
        LOG << LOGMSG << "collate starts." << endl;
        mr2->collate(NULL);
        LOG << LOGMSG << "collate ends." 
            << "\t" << MPI_Wtime() - collate_time << endl;
        
        double reduce_time = MPI_Wtime();
        LOG << LOGMSG << "reduce starts." << endl;
        mr2->reduce(&mr_sort_multivalues_by_evalue, NULL);
        LOG << LOGMSG << "reduce ends." 
            << "\t" <<  MPI_Wtime() - reduce_time << endl;
        
        ///
        /// Save history
        ///
        if (MYID == 0) {                
            string histFileName = OUTPREFIX + "-history.txt";
            ofstream histFile(histFileName.c_str(), ios::out | ios::app);
            if (!histFile) {
                cerr << "ERROR: failed to open a history file" << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else {
                time_t timer;
                timer = time(NULL);
                histFile << subMasterFileName << "," 
                         << asctime(localtime(&timer));
            }
            histFile.close();
        }
            
        MPI_Barrier(MPI_COMM_WORLD);           
    }
     
    delete mr2;
    delete pTargetDb;      
    if (LOGORNOT && LOGTOFILE) logFile.close();    
    
    profile_time = MPI_Wtime() - profile_time;
    if (MYID == 0) cout << "Total Execution Time: " << profile_time << endl;
    
    MPI_Finalize();
            
    return 0;
}

/** MR-MPI Map function - settup NCBI C++ Toolkit env and call blast
 * @param itask
 * @param file
 * @param kv
 * @param ptr
 */

void mr_run_blast(int itask, 
                  char *file, 
                  KeyValue *kv, 
                  void *ptr)
{
    /// 
    /// Make a option handle and cblastinputsource
    ///
    //CRef<CBlastnAppArgs> CmdLineArgs; 
    //CmdLineArgs.Reset(new CBlastnAppArgs());
    
    /*** Get the BLAST options ***/
    //const CArgs& args = GetArgs();
    //string allArgs;
    //cout << "All args = " << args.Print(allArgs) << endl;
    //import_search_strategy(STRATEGYFILENAME, CmdLineArgs);
    
    ifstream strategyFile(STRATEGYFILENAME.c_str(), ios::in);
    if (!strategyFile.is_open()) {
        cerr << "ERROR: failed to open a search strategy file" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    CNcbiIstream* in = &strategyFile;
    CRef<CBlast4_request> b4req;
    try { 
        b4req = ExtractBlast4Request(*in);
    } catch (const CSerialException&) {
        NCBI_THROW(CInputException, eInvalidInput, 
                   "Failed to read search strategy");
    }
    CImportStrategy strat(b4req);
    CRef<blast::CBlastOptionsHandle> opts = strat.GetOptionsHandle();
    
    //CRef<CBlastOptionsHandle> opts_hndl(&*CmdLineArgs->SetOptions(args));
    const CBlastOptions& opt = opts->GetOptions();
    
    //EProgram program = ProgramNameToEnum("blastn");
    //CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    //set_default_opts(opts);    
    //opts->Validate();

    if (MYID == 1 && OPTDUMPED == false) {
        OPTDUMPED = true;
        opts->GetOptions().DebugDumpText(LOG, "opts", 1);
    }

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    
    ///
    /// Split work item: <beginOffset,endOffset,db chunk name>
    ///
    string workItem(file);
    vector<string> vWorkItemTokens;
    boost::split(vWorkItemTokens, workItem, boost::is_any_of(","));
    assert(vWorkItemTokens.size() == 3);

    ///
    /// Target db chunk setting
    ///
    string dbChunkName = vWorkItemTokens[2];
    LOG << LOGMSG << "DB name = " << dbChunkName << endl;
    
    double db_load_time = MPI_Wtime();
    LOG << LOGMSG << "DB loading starts." << endl;
    if(pTargetDb == 0 || dbChunkName != prevDbChunkName) {
        delete pTargetDb;
        pTargetDb = new CSearchDatabase(dbChunkName,
                                        CSearchDatabase::eBlastDbIsNucleotide);
    }
    prevDbChunkName = dbChunkName;
    LOG << LOGMSG << "DB loading ends." 
        << "\t" << MPI_Wtime() - db_load_time << endl;

    ///
    /// Read sequence block and run blast
    ///
    unsigned long long beginOffset 
        = boost::lexical_cast<unsigned long long>(vWorkItemTokens[0]);
    unsigned long long endOffset 
        = boost::lexical_cast<unsigned long long>(vWorkItemTokens[1]);
    
    LOG << LOGMSG << "Offsets = " << beginOffset << " " << endOffset;
    
    unsigned long long blSize = endOffset-beginOffset;
    char* buff2 = (char*) malloc(sizeof(char)*blSize);
    FILE *qf = fopen((char*)QUERYFILENAME.c_str(), "r");
    
    fseek(qf, beginOffset, SEEK_SET);    
    string query;
    vector<string> vHeaders;
    unsigned int nQuery = 0;
    while (ftell(qf) < endOffset) {
        fgets(buff2, blSize, qf);
        if (buff2[0] == '>') {
            vHeaders.push_back(string(buff2));
            LOG << " " << string(buff2) << endl;
            nQuery++;
        }
        query += string(buff2);
    }    
    fclose(qf);
    free(buff2);
    LOG << LOGMSG << "Num query for Blast call = " << nQuery << endl;        
 
    ///
    /// Set queries as fasta input
    ///
    CBlastFastaInputSource fasta_input(query, iconfig);
    CBlastInput blastInput(&fasta_input);
    CScope scope(*objmgr);
    TSeqLocVector queryLoc = blastInput.GetAllSeqLocs(scope);
    CRef<IQueryFactory> queryFactory(new CObjMgr_QueryFactory(queryLoc));

    ///
    /// Run blast
    ///
    double blaster_init_time = MPI_Wtime();
    LOG << LOGMSG << "Blaster init starts." << endl;
    CLocalBlast blaster(queryFactory, opts, *pTargetDb);
    LOG << LOGMSG << "Blaster init ends." 
           << "\t" << MPI_Wtime() - blaster_init_time << endl;
        
    double blast_call_time = MPI_Wtime();
    
    LOG << LOGMSG << "Blast call starts." << endl;
    //////////////////////////////////////////
    CSearchResultSet results = *blaster.Run();
    //////////////////////////////////////////
    LOG << LOGMSG << "Blast call ends." 
           << "\t" << MPI_Wtime() - blast_call_time << endl;
    
    ///
    /// Get warning messages
    ///
    for (uint64_t i = 0; i < results.GetNumResults(); ++i) {
        TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
        if (messages.size() > 0) {
            CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            if (seq_id.NotEmpty())
                cerr << "ID: " << seq_id->AsFastaString() << endl;
            else
                cerr << "ID: " << "Unknown" << endl;

            ITERATE(vector<CRef<CSearchMessage> >, it, messages) {
                cerr << (*it)->GetMessage() << endl;
            }
        }
    }

    ///
    /// Get the results
    ///   
    double adding_kv_time = MPI_Wtime();
    LOG << LOGMSG << "Adding hits to KV starts." << endl;
    
    for (uint64_t i = 0; i < results.GetNumResults(); ++i) {

        CConstRef<CSeq_align_set> aln_set = results[i].GetSeqAlign();

        if (results[i].HasAlignments()) {
            ITERATE(CSeq_align_set::Tdata, itr, aln_set->Get()) {
                const CSeq_align& s = **itr;
                
                ///
                /// VERY IMPORTANT!!
                /// queryID is not unique. It's internal qid in fasta_input.
                ///
                string queryID = s.GetSeq_id(QUERY).GetSeqIdString(); 
                string subID = s.GetSeq_id(SUBJECT).GetSeqIdString(); /// GI
                
                ///
                /// Refer: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/
                /// doxyhtml/classCSeq__align.html
                ///
                double pIdentity=0.0;
                s.GetNamedScore(CSeq_align::eScore_PercentIdentity, 
                    pIdentity);

                int alignLen=0, qStart=0, qEnd=0;
                int sStart=0, sEnd=0;
                int gapOpens=0;
                gapOpens    = s.GetNumGapOpenings();
                qStart      = s.GetSeqStart(QUERY);
                qEnd        = s.GetSeqStop(QUERY);
                sStart      = s.GetSeqStart(SUBJECT);
                sEnd        = s.GetSeqStop(SUBJECT);
                alignLen    = s.GetAlignLength();
                
                double eValue=0.0;
                int bitScore=0;
                int misMatches=0;
                
                s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);
                s.GetNamedScore(CSeq_align::eScore_MismatchCount, misMatches);
                
                ///
                /// Tokenize query header
                /// 
                string qHeader = 
                    vHeaders[boost::lexical_cast<unsigned int>(queryID)-1];
                vector<string> vQueryId;
                boost::split(vQueryId, qHeader, boost::is_any_of("|"));
                string uniqueQID = vQueryId[2]; /// query id
                
                /// 
                /// Add a csv blast result to kv
                /// 
                if (EXCLUSIONORNOT) {     
                    ///
                    /// The below infomation is only for our own simulated
                    /// sequence data sets which has 
                    /// - origin GI
                    /// - unique query ID                        
                    /// - original sequence length
                    /// - cut location start
                    /// - cut location end
                    ///
                    string qGi = vQueryId[1];   /// GI
                    //uint64_t origLen            /// length of the orig seq         
                        //= boost::lexical_cast<unsigned int>(vQueryId[3]);     
                    int qCutLocStart            /// cut coordinates - start
                        = boost::lexical_cast<unsigned int>(vQueryId[4]);      
                    uint64_t qCutLocEnd         /// cut coordinates - end
                        = boost::lexical_cast<unsigned int>(vQueryId[5]);  
                                           
                    if (!check_exclusion(qGi, subID, qCutLocStart, qCutLocEnd, 
                        sStart, sEnd, EXCLUSIONTHRESHOLD)) {
                                                
                        ///
                        /// To pass Blast hits using struct, outfmt=6
                        /// query id, subject id, % identity, alignment length, 
                        /// mismatches, gap opens, q. start, q. end, s. start, 
                        /// s. end, evalue, bit score
                        ///
                        BLASTRES res;
                        res.subjectid 
                            = boost::lexical_cast<unsigned long long>(subID);
                        res.identity = pIdentity;
                        res.alignlen = alignLen;
                        res.mismatches = misMatches;
                        res.gapopens = gapOpens;
                        res.qstart = qStart;
                        res.qend = qEnd;
                        res.sstart = sStart;
                        res.send = sEnd;
                        res.evalue = eValue;
                        res.bitscore = bitScore;     
                                                
                        /// 
                        /// ADD <KEY = "QUERYID", VALUE="BLASTRESULT">
                        /// TO KV
                        ///
                        const char* newKey = (char*)((uniqueQID).c_str()); 
                        kv->add((char*)newKey, strlen(newKey) + 1, 
                            (char*)&res, BLASTRESSZ);
                    
                    }
                    /// 
                    /// Found a self hits. Record the hits in a histroy file
                    ///
                    else {
                        string exFileName = OUTPREFIX + "-" 
                            + boost::lexical_cast<string>(MYID) + "-" 
                            + EXCLUSIONHISTFILE;
                        ofstream exFile(exFileName.c_str(), ios::out 
                            | ios::app);
                            
                        if (!exFile) {
                            cerr << "ERROR: failed to open a exclusion "
                                 << "history file" << endl;
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                        else {
                            /// Format:
                            /// qGi sGi qCutLocStart qCutLocEnd sStart sEnd
                            exFile << qGi << "\t"
                                   << subID << "\t"
                                   << qCutLocStart << "\t"
                                   << qCutLocEnd << "\t" 
                                   << sStart << "\t"
                                   << sEnd << endl;
                        }
                        exFile.close();
                    }
                }
                else {
                    BLASTRES res;
                    res.subjectid 
                        = boost::lexical_cast<unsigned long long>(subID);
                    res.identity    = pIdentity;
                    res.alignlen    = alignLen;
                    res.mismatches  = misMatches;
                    res.gapopens    = gapOpens;
                    res.qstart      = qStart;
                    res.qend        = qEnd;
                    res.sstart      = sStart;
                    res.send        = sEnd;
                    res.evalue      = eValue;
                    res.bitscore    = bitScore;     
                                            
                    const char* newKey = (char*)((uniqueQID).c_str()); 
                    kv->add((char*)newKey, strlen(newKey) + 1, (char*)&res, 
                        BLASTRESSZ);
                }
            }
        }
    }
    LOG << LOGMSG << "Adding hits to KV ends." 
           << "\t" << MPI_Wtime() - adding_kv_time << endl;
}


/** Sort function - Passed to MR-MPI sort_values() for sorting blast result
 * string by bit score.
 * @param e1
 * @param e2
 */
 
inline bool mycompare(STRUCTEVALUE e1, 
                      STRUCTEVALUE e2)
{
    return (e1.evalue < e2.evalue);
}

/** Sort by wvalue - Passed to MR-MPI reduce() for sorting KMVs by evalue.
 * @param key
 * @param keybytes
 * @param multivalue - collected blast result strings. There could be more than
 * nvalues results in it. The separator '>' should be used for splitting the 
 * resutls.
 * @param nvalues
 * @param valuebytes
 * @param kv
 * @param ptr
 */
 
void mr_sort_multivalues_by_evalue(char *key, 
                                   int keybytes, 
                                   char *multivalue,
                                   int nvalues, 
                                   int *valuebytes, 
                                   KeyValue *kv, 
                                   void *ptr) 
{
    double sort_and_save_time = MPI_Wtime();
    LOG << LOGMSG << "Sort/save starts." << endl;
    
    /// Check if there is KMV overflow
    assert(multivalue != NULL && nvalues != 0);
    
    ///
    /// Make STRUCTEVALUE = {BLASTRES* p; double evalue;}
    /// and sort by evalue
    ///
    vector<STRUCTEVALUE> vforsort;
    for (int n = 0; n < nvalues; n++) {                
        BLASTRES* res = (BLASTRES*)multivalue;
        STRUCTEVALUE structEvalue;
        structEvalue.p = res;
        structEvalue.evalue = res->evalue;
        vforsort.push_back(structEvalue);
        multivalue += BLASTRESSZ;         
    }
    sort(vforsort.begin(), vforsort.end(), mycompare);
    
    ///
    /// After sorting, each worker saves the set of results into a file.
    /// Note: The file open option is "a".
    /// Note: outfmt = 6 in Blast
    ///
    ofstream outputFile(WORKEROUTPUTFILENAME.c_str(), ios::out | ios::app);
    for (int n = 0; n < nvalues; n++) {
        BLASTRES* res = (BLASTRES*)(vforsort[n].p);
        outputFile << key << "\t"
            << res->subjectid << "\t"
            << res->identity << "\t"
            << res->alignlen << "\t"
            << res->mismatches << "\t"
            << res->gapopens << "\t"
            << res->qstart << "\t"
            << res->qend << "\t"
            << res->sstart << "\t"
            << res->send << "\t"
            << res->evalue << "\t"
            << res->bitscore 
            << endl;
    }
    outputFile.close();
    vforsort.clear();
    
    LOG << LOGMSG << "Sort/save ends." 
           << "\t" << MPI_Wtime() - sort_and_save_time << endl;
}   
 
/** Set Blast options
 * @param optsHandle
 */
 
void set_default_opts(CRef<CBlastOptionsHandle> optsHandle)
{
    optsHandle->SetDefaults();
    optsHandle->SetEvalueThreshold(EVALUE);
    //optsHandle->SetMatchReward(0);
    //optsHandle->SetMismatchPenalty(0);
    //optsHandle->SetMatrixName("BLOSUM62");
    //optsHandle->SetCutoffScore(CUTOFFSCORE);
    if (CBlastNucleotideOptionsHandle* nuclHandle =
                dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {

        //nuclHandle->SetMatchReward(0);
        //nuclHandle->SetMismatchPenalty(0);
        nuclHandle->SetEvalueThreshold(EVALUE);
    }

    return;
}


void import_search_strategy(string strategyFile, 
                            blast::CBlastAppArgs* cmdline_args)
{
    ////CNcbiIstream* in = cmdline_args->GetImportSearchStrategyStream(args);
    ifstream s(strategyFile.c_str(), ios::in);
    CNcbiIstream* in = &s;
    if ( !in ) {
        return;
    }
    const bool is_remote_search = false;
    const bool override_query = true;
    const bool override_subject = true;
    s_ImportSearchStrategy(in, cmdline_args, is_remote_search, override_query,
                           override_subject);
    //if (CMbIndexArgs::HasBeenSet(args)) {
        //ERR_POST(Warning << "Overriding megablast BLAST DB indexed options in saved strategy");
    //}
}

/// Imports search strategy, using CImportStrategy.
static void s_ImportSearchStrategy(CNcbiIstream* in, 
                                   blast::CBlastAppArgs* cmdline_args,
                                   bool is_remote_search, 
                                   bool override_query, 
                                   bool override_subject)
{
    //if ( !in ) {
        //return;
    //}

    //CRef<CBlast4_request> b4req;
    //try { 
        //b4req = ExtractBlast4Request(*in);
    //} catch (const CSerialException&) {
        //NCBI_THROW(CInputException, eInvalidInput, 
                   //"Failed to read search strategy");
    //}

    //CImportStrategy strategy(b4req);

    //CRef<blast::CBlastOptionsHandle> opts_hndl = strategy.GetOptionsHandle();
    //cmdline_args->SetOptionsHandle(opts_hndl);
    //const EBlastProgramType prog = opts_hndl->GetOptions().GetProgramType();
    //cmdline_args->SetTask(strategy.GetTask());

    //// Get the subject
    //if (override_subject) {
        //ERR_POST(Warning << "Overriding database/subject in saved strategy");
    //} else {
        //CRef<blast::CBlastDatabaseArgs> db_args;
        //CRef<CBlast4_subject> subj = strategy.GetSubject();
    //const bool subject_is_protein = Blast_SubjectIsProtein(prog) ? true : false;

        //if (subj->IsDatabase()) {
            //CBlastOptionsBuilder bob(strategy.GetProgram(), strategy.GetService(), CBlastOptions::eBoth);
            //bob.GetSearchOptions(&strategy.GetAlgoOptions(), &strategy.GetProgramOptions());
            //db_args = s_ImportDatabase(*subj, bob, subject_is_protein,
                                       //is_remote_search);
        //} else {
            //db_args = s_ImportSubjects(*subj, subject_is_protein);
        //}
        //_ASSERT(db_args.NotEmpty());
        //cmdline_args->SetBlastDatabaseArgs(db_args);
    //}

    //// Get the query, queries, or pssm
    //if (override_query) {
        //ERR_POST(Warning << "Overriding query in saved strategy");
    //} else {
        //CRef<CBlast4_queries> queries = strategy.GetQueries();
        //if (queries->IsPssm()) {
            //s_ImportPssm(*queries, opts_hndl, cmdline_args);
        //} else {
            //s_ImportQueries(*queries, opts_hndl, cmdline_args);
        //}
        //// Set the range restriction for the query, if applicable
        //const TSeqRange query_range = strategy.GetQueryRange();
        //if (query_range != TSeqRange::GetEmpty()) {
            //cmdline_args->GetQueryOptionsArgs()->SetRange(query_range);
        //}
    //}
}

/** Check exclusion - Based on the coordinates of query and subject, decide 
 * whether the result should be included in the final result or not.
 * @param qGi: query ID
 * @param sGi: subject ID
 * @param qCutLocStart: Cutting start location of the query from original 
 * refseq_genomic fasta input.
 * @param qCutLocEnd: Cutting end location of the query from original 
 * refseq_genomic fasta input.
 * @param sStart: subject alignment start location.
 * @param sEnd: subject alignemnt end location.
 * @param threshold: overlap threshold (default = 100bp).
 */
 
bool check_exclusion(string qGi, 
                     string sGi, 
                     uint64_t qCutLocStart,
                     uint64_t qCutLocEnd, 
                     uint64_t sStart, 
                     uint64_t sEnd,
                     int threshold)
{
    /// 
    /// To exclude Blast result from the original sequence from which
    /// the input query is originated (sampled). Basically if qGi == sGi
    /// and qCutLocStart is similar with sStart and qCutLocEnd is similar
    /// with sEnd in terms of coordinates, the result should be excluded.
    /// 
    /// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    ///                           |              |
    ///                     qCutLocStart      qCutLocEnd
    /// 
    /// Query:                    XXXXXXXXXXXXXXXX
    ///                              |          |
    ///                           qStart       qEnd
    /// 
    /// Subject:   ------------------XXXXXXXXXXXX-----------------------
    ///                              |          |
    ///                           sStart      sEnd
    /// 
    bool ret = false;

    if (qGi == sGi) {
        if (qCutLocStart < 0) {
            /// 
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            /// 
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }
        if ((qCutLocStart - threshold <= sStart 
            && sStart <= qCutLocStart + threshold) 
            && (qCutLocEnd - threshold <= sEnd 
            && sEnd <= qCutLocEnd + threshold))
            ret = true;
    }

    return ret;
}

/// EOF

