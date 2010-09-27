////////////////////////////////////////////////////////////////////////////////
//
//  Parallelizing BLAST on MR-MPI
//
//  Author: Seung-Jin Sul
//          (ssul@jcvi.org)
//
//  Revisions
//
//      09.01.2010  Start!
//
//      09.05.2010  blast call done
//
//      09.07.2010  Adding KV in a mapper and aggregating the vecResults
//
//      09.09.2010  Start v2. Decide to run several MR-MPI jobs and each
//                  job searches on a specific DB chunk.
//
//      09.10.2010  v2 done!
//                  Add print2file member function to KeyValue class.
//
//      09.15.2010  v3 with NCBI C++ Toolkit started.
//
//      09.20.2010  Running version is done. ASN.1 parsing should be done.
//
//      09.21.2010  XML parsing tried but returned back to use raw blast 
//                  output for creating KV
//
//      09.22.2010  Using raw blast output, KV is constructed. Running 
//                  version is done. Start testing on Ranger.
//
////////////////////////////////////////////////////////////////////////////////

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

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
#include "./mrmpi/mapreduce.h"
#include "./mrmpi/keyvalue.h"

using namespace MAPREDUCE_NS;
using namespace std;

const size_t MAXSTR         = 256;
const size_t MAXQUERYNUM    = 1000;

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

/// ASN
#include <sstream>

/// blastn importing
//#include <algo/blast/blastinput/blastn_args.hpp>
//#include <algo/blast/format/blast_format.hpp>
//#include "blast_app_util.hpp"

/// Tabular
#include <algo/blast/format/blast_format.hpp>
#include <algo/blast/format/blastxml_format.hpp>
#include <algo/blast/format/data4xmlformat.hpp>       /* NCBI_FAKE_WARNING */
#include <algo/blast/format/build_archive.hpp>
#include <objects/seq/Seq_annot.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>
#include <algo/blast/core/blast_stat.h>
#include <corelib/ncbiutil.hpp>                 // for FindBestChoice

/// XML
//#include <misc/xmlwrapp/xmlwrapp.hpp>           // for XML::
//#include <misc/xmlwrapp/xsltwrapp.hpp>

/// CSeq_align
#include <objects/seqalign/Seq_align.hpp>
#include <util/range.hpp>


USING_NCBI_SCOPE;
USING_SCOPE(blast);
//USING_SCOPE(objects);
//USING_SCOPE(align_format);
//USING_SCOPE(xml);

const size_t QUERY          = 0;
const size_t SUBJECT        = 1;


/// TO COLLECT PROC NAMES ON WHICH EACH TASK IS ASIGNED
//struct procName {
    //int     myId;
    //char*   pName;
//};


/// FOR COMMUNICATING WITH MR-MPI
struct giftBox {
    int             myId;
    int             numProc;
    //int             numDbChunks;
    //int             totalDbChunks;
    char*           pName;
    //char*           queryFileName;
    char*           dbChunkName;
    //unsigned long   numQuery;
    //CBlastDemoApplication* cbDemo;
    //multimap<string, unsigned int> *mmapPName;
    //vector<string>* vPName;
};

/// MR-MPI CALLS
//void    get_node_name(int, KeyValue *, void *);
//void    collect_node_names(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
//void    save_node_names(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void    set_default_opts(CRef<CBlastOptionsHandle> optsHandle);
void    run_blast(int itask, char *file, KeyValue *kv, void *ptr);
//int     key_compare_str(char *p1, int len1, char *p2, int len2);
int     key_compare_int(char *p1, int len1, char *p2, int len2);

bool    check_exclusion(string qGi, string sGi, int qCutLocStart, size_t qCutLocEnd, size_t sStart, size_t sEnd, size_t threshold);
    
/// TOKENIZER ROUTINES
vector<string> &split(const string &s, char delim, vector<string> &vecElems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        vecElems.push_back(item);
    }
    return vecElems;
}

vector<string> split(const string &s, char delim)
{
    vector<string> vecElems;
    return split(s, delim, vecElems);
}

/// NOTE: NCBI PROVIDES THESE UTILS
/// INT TO STRING 
string int2str(int number)
{
    stringstream ss;
    ss << number;
    return ss.str();
}

size_t str2uint(string str)
{
    std::stringstream ss;
    ss << str;
    size_t f;
    ss >> f;
    return f;
}

int str2int(string str)
{
    std::stringstream ss;
    ss << str;
    int f;
    ss >> f;
    return f;
}



int main(int argc, char* argv[])
{
    // Execute main application function
    //int newArgc = 7;
    //char *newArgv[] = {"./mrblast", "-program", "blastn", "-db", 
        //"nt.00", "-in", "test2.query"};
    //CBlastDemoApplication demo;
    //demo.AppMain(newArgc, newArgv);

    /// MPI
    MPI_Init(&argc, &argv);

    char MPI_procName[100];

    int MPI_myId, MPI_nProcs, MPI_length;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myId);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProcs);
    MPI_Get_processor_name(MPI_procName, &MPI_length);
    fprintf(stdout, "### [Node %d]: %s \n", MPI_myId, MPI_procName);

    //if (argc != 4) {
    //if (MPI_myId == 0) printf("Syntax: mpirun -np n mrblast masterFileName 
        //chunkName outFileName\n");
    //MPI_Abort(MPI_COMM_WORLD, 1);
    //}

    /*
    ///MR-MPI for node name collecting
    MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
    mr->verbosity = 0;
    mr->timer = 1;

    MPI_Barrier(MPI_COMM_WORLD);
    //double tstart = MPI_Wtime();

    ///
    /// GET MPI PROC NAMES
    /// TO ASSIGN A DB CHUNK WITH A NODE SO THAT THE TASKS ON THE NODE
    /// SEARCHES QUERIES ON THE DB CHUNK
    ///
    procName nn;
    nn.myId = MPI_myId;
    nn.pName = MPI_procName;
    int sizeKV = mr->map(MPI_nProcs, &get_node_name, &nn);
    //mr->print(-1, 1, 5, 5);

    /// DEBUG
    //mr->collate(NULL);
    //mr->print(-1, 1, 5, 5);
    //int numUniqKV = mr->reduce(&collect_node_names, (void *)NULL);
    //mr->print(-1, 1, 5, 5);

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    mr->gather(1); /// KVs FROM TASKs -> A KV IN THE ROOT NODE
    mr->print(-1, 1, 5, 5);

    mr->broadcast(0);
    //mr->print(-1, 1, 5, 5);

    //multimap<string, unsigned int> mmapProcNames;
    vector<string> vecProcName;
    mr->map(mr, &save_node_names, &vecProcName);
    //cout << vecProcName.size() << endl;

    MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    /// DEBUG
    //if (MPI_myId == 0) {
    ////pair<multimap<string, unsigned int>::iterator, multimap<string, unsigned int>::iterator> ppp;
    ////ppp = mmapProcNames.equal_range("ssul-lx");
    ////for (multimap<string, unsigned int>::iterator it2 = ppp.first; it2 != ppp.second; ++it2) {
    ////cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << endl;
    ////}
    //for (size_t i = 0; i < vecProcName.size(); ++i)
    //cout << "### Node " << MPI_myId << " name = " << vecProcName[i] << endl;
    //}

    delete mr;
    */
    
    ///
    /// Blast search
    ///
    MapReduce *mr2 = new MapReduce(MPI_COMM_WORLD);
    mr2->verbosity = 0;
    mr2->timer = 1;
    mr2->mapstyle = 2; /// master/slave mode

    int nvecRes;
    char *masterFileName = argv[1];

    giftBox gf;
    gf.myId = MPI_myId;
    //gf.numProc = MPI_nProcs;
    gf.pName = MPI_procName;
    gf.dbChunkName = argv[2];
    //gf.mmapPName = &mmapProcNames;
    //gf.vPName = &vecProcName;
    //gf.cbDemo = &demo;
    
    //nvecRes = mr2->map(masterFileName, &run_blast, &gf);
    //cout << "### [Node " << MPI_myId << "] nvecRes = " << nvecRes << endl;
    //mr2->print(-1, 1, 5, 5);
    
    //mr2->collate(NULL); //////////////////////////////////////////////
    //mr2->print(-1, 1, 5, 5);

    //MPI_Barrier(MPI_COMM_WORLD); ///////////////////////////////////////

    //mr2->gather(1);
    //mr2->print(0, 1, 5, 5);

    //mr2->sort_keys(&key_compare_int);
    //mr2->print(0, 1, 5, 5);
    
    //if (MPI_myId == 0) {
        //cout << "Saving...\n";
        ////nvecRes = mr2->map(mr2, &output, &gf);
        //string outFileName = string(argv[3]) + ".txt";
        //FILE *outFile = fopen(outFileName.c_str(), "w");
        //if (outFile) {
            //mr2->kv->print2file(1, 5, 5, outFile);
        //}
        //cout << "DONE!\n";
    //}

    delete mr2;
    MPI::Finalize();
    cout << "Finalized\n";
}


/* ---------------------------------------------------------------------- */
void run_blast(int itask, char *file, KeyValue *kv, void *ptr)
/* ---------------------------------------------------------------------- */
{
    giftBox *gf = (giftBox *) ptr;

    ///
    /// Make args, argv for AppMain
    /// ex) blast -program blastn -db nt.00 -in test.query
    ///
    //int newArgc = 7;
    //char *newArgv[] = {"./mrblast", "-program", "blastn", "-db", 
        //"nt.00", "-in", "test2.query"};
    //CBlastDemoApplication().AppMain(newArgc, newArgv);
    //demo->AppMain(newArgc, newArgv);

    ///
    /// Make a option handle and CBlastInputSource
    ///
    EProgram program = ProgramNameToEnum("blastn");
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    set_default_opts(opts);
    opts->Validate();

    /// DEBUG
    //if (gf->myId == 0) 
        //opts->GetOptions().DebugDumpText(cerr, "opts", 1);

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        throw std::runtime_error("Could not initialize object manager");
    }

    bool isProtein = false;
    SDataLoaderConfig dlconfig(isProtein);
    CBlastInputSourceConfig iconfig(dlconfig);
    ifstream queryFile(file);
    //ifstream queryFile("test2.query");
    cout << "qFile Name = " << file << endl;
    //string qFileName = "test2.query";
    //ifstream queryFile(qFileName.c_str());
    //CBlastFastaInputSource fasta_input(qFileName, iconfig);

    /// Use query string instead of file input
    string header, seq;
    vector<string> vecSeq;

    ///
    /// Target DB setting
    ///
    cout << "dbChunkName = " << gf->dbChunkName << endl;
    string dbChunkName(gf->dbChunkName);
    const CSearchDatabase target_db(dbChunkName, CSearchDatabase::eBlastDbIsNucleotide);
    
    ///
    /// Read seq(s) and run blast
    ///
    while (!getline(queryFile, header).eof()) {
        getline(queryFile, seq);
        
        vector<string> vHeaders;
        
        if (header.length() == 0 || seq.length() == 0) {
            fprintf(stderr, "### ERROR: [Node %d]: %s, itask = %d, \
                file = %s, Seq read error.\n",
                gf->myId, gf->pName, itask, file);
            exit(0);
        }
        fprintf(stdout, "### INFO: [Node %d]: %s, itask = %d, %s\n",
            gf->myId, gf->pName, itask, header.c_str());

        string query = header + '\n' + seq;
        vHeaders.push_back(header);
        
        ///
        /// Accummulate seqs when MAXQUERYNUM > 1
        ///
        size_t j = 1;
        while (j < MAXQUERYNUM && !getline(queryFile, header).eof()) {
            getline(queryFile, seq);
            fprintf(stdout, "### INFO: [Node %d]: %s, itask = %d, %s, num queries = %d\n",
                gf->myId, gf->pName, itask, header.c_str(), j + 1);
            query = query + '\n' + header + '\n' + seq;
            vHeaders.push_back(header);
            j++;
        }
        
        /// Tokenize query header
        //vector<string> vQeuryId  = split(header, '|');
        //string gi                = vQeuryId[1];             /// GI
        //string qid               = vQeuryId[2];             /// QUERY ID
        //size_t oLen              = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
        //int qStart               = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
        //size_t qEnd              = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
        //fprintf(stderr, "### INFO: [Node %d]: %s, itask = %d, %s %s %d %d %d\n",
                //gf->myId, gf->pName, itask, gi.c_str(), qid.c_str(), oLen, qStart, qEnd);
                        
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
        CLocalBlast blaster(queryFactory, opts, target_db);
        CSearchResultSet results = *blaster.Run();

        ///
        /// Get the results
        ///

        /// Get warning messages.
        for (size_t i = 0; i < results.GetNumResults(); i++) {
            TQueryMessages messages = results[i].GetErrors(eBlastSevWarning);
            if (messages.size() > 0) {
                CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
                if (seq_id.NotEmpty())
                    cerr << "ID: " << seq_id->AsFastaString() << endl;
                else
                    cerr << "ID: " << "Unknown" << endl;

                ITERATE(vector<CRef<CSearchMessage> >, it, messages)
                cerr << (*it)->GetMessage() << endl;
            }
        }

        cout << "results.GetNumResults() = " << results.GetNumResults() << endl;
        cout << "results.GetNumQueries() = " << results.GetNumQueries() << endl;
        for (size_t i = 0; i < results.GetNumResults(); i++) {

            //CConstRef<CSeq_id> seq_id = results[i].GetSeqId();
            //cout << "qID = " << seq_id->AsFastaString() << endl;

            CConstRef<CSeq_align_set> aln_set = results[i].GetSeqAlign();
            //cout << MSerial_AsnText << *aln_set;
            //cout << MSerial_Xml << *aln_set;

            /// XML parsing
            /*
            stringstream oss;
            oss << MSerial_Xml << *aln_set;
            //std::cerr << "[" << oss.str() << "]\n";
            string xmldata(oss.str());
            cout << xmldata<< endl;
            
            cout << "xmldata.size() = " << xmldata.size() << endl;
            //std::string xmldata( "<TagA>"
                    //"<TagB>stuff</TagB>"
                    //"</TagA>" );
            xml::tree_parser parser(xmldata.c_str(), xmldata.size());
            //xml::tree_parser parser( "xml-ex.txt" );
            xml::document &doc = parser.get_document();
            const xml::attributes & attrs =
                parser.get_document().get_root_node().get_attributes();
            xml::attributes::const_iterator i = attrs.begin();
            xml::attributes::const_iterator end = attrs.end();
            for (; i != end; ++i) {
                std::cout << i->get_name() << "=" << i->get_value() << "\n";
            }
            
            xml::node & root = parser.get_document().get_root_node();
            xml::node::const_iterator child( root.begin() );
            xml::node::const_iterator child_end( root.end() );
            std::cout << "root node is '" << root.get_name() << "'\n";
            for ( ; child != child_end; ++child ) {
                if ( child->is_text() ) continue;
                std::cout << "child node '" << child->get_name() << "'\n"; 
            }
            //xml::ns rootSpace( root.get_namespace() );
            //std::cout << "Root namespace: " << rootSpace.get_prefix() << "->"
                //<< rootSpace.get_uri() << "\n";
            */
            
            ///*
            /// Print tabular format
            //CNcbiOstream& out = args["out"].AsOutputFile();
            //CBioseq_Handle bhandle = scope.GetBioseqHandle(results[i].GetSeqId(),
                                                      //CScope::eGetBioseq_All);

            //const CBlastTabularInfo::EFieldDelimiter kDelim = CBlastTabularInfo::eComma;
            ////(m_FormatType == CFormattingArgs::eCommaSeparatedValues
             ////? CBlastTabularInfo::eComma : CBlastTabularInfo::eTab);
            //CBlastTabularInfo tabinfo(cout, "", kDelim);
            //tabinfo.SetParseLocalIds(true);
            //CNcbiMatrix<int> m_ScoringMatrix;
            //CAlignFormatUtil::GetAsciiProteinMatrix("BLOSUM62", m_ScoringMatrix);

            //CConstRef<CBioseq> subject_bioseq = x_CreateSubjectBioseq();
            //tabinfo.PrintHeader(strProgVersion, *(bhandle.GetBioseqCore()),
                                //m_DbName, results.GetRID(), itr_num, aln_set,
                                //subject_bioseq);
                                            
            if (results[i].HasAlignments()) {
                ITERATE(CSeq_align_set::Tdata, itr, aln_set->Get()) {
                    const CSeq_align& s = **itr;
                    //tabinfo.SetFields(s, scope, &m_ScoringMatrix);
                    //tabinfo.Print();
                    
                    string qNum = s.GetSeq_id(QUERY).GetSeqIdString();           
                    //cout << s.GetSeq_id(QUERY).AsFastaString() << endl;
                    //cout << s.GetSeq_id(QUERY).GetSeqIdString() << endl;
                    string seqId = s.GetSeq_id(SUBJECT).GetSeqIdString();
                    
                    //double pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage;
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Gapped, pIdentityGapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_Ungapped, pIdentityUngapped);
                    //s.GetNamedScore(CSeq_align::eScore_PercentIdentity_GapOpeningOnly, pIdentityGapOpeningOnly);
                    //s.GetNamedScore(CSeq_align::eScore_PercentCoverage, pCoverage);
                    
                    size_t alignLen, qStart, qEnd, sStart, sEnd;
                    //CRange<TSeqPos> qRange, sRange;
                    //qRange = s.GetSeqRange(QUERY);
                    //sRange = s.GetSeqRange(SUBJECT);
                    qStart      = s.GetSeqStart(QUERY);
                    qEnd        = s.GetSeqStop(QUERY);
                    sStart      = s.GetSeqStart(SUBJECT);
                    sEnd        = s.GetSeqStop(SUBJECT);
                    alignLen    = s.GetAlignLength();
                    
                    double eValue;
                    int genericScore, bitScore;
                    s.GetNamedScore(CSeq_align::eScore_EValue, eValue);
                    s.GetNamedScore(CSeq_align::eScore_Score, genericScore);                    
                    s.GetNamedScore(CSeq_align::eScore_BitScore, bitScore);       
                    
                    ///
                    /// Tokenize query header
                    ///
                    string qHeader = vHeaders[str2uint(qNum)-1];
                    //cout << "qHeader = " << qHeader << endl;
                    vector<string> vQeuryId  = split(qHeader, '|');
                    string qGi               = vQeuryId[1];             /// GI
                    string qid               = vQeuryId[2];             /// QUERY ID
                    size_t origLen           = str2uint(vQeuryId[3]);   /// LENGTH OF THE ORIG SEQ
                    int qCutLocStart         = str2int(vQeuryId[4]);    /// CUT COORDINATES - START
                    size_t qCutLocEnd        = str2uint(vQeuryId[5]);   /// CUT COORDINATES - END
                    fprintf(stderr, "### INFO: [Node %d]: %s, itask = %d, %s %s %d %d %d\n",
                        gf->myId, gf->pName, itask, qGi.c_str(), qid.c_str(), 
                        origLen, qCutLocStart, qCutLocEnd);
                        
                    ///
                    /// Add a CSV blast result to KV
                    ///
                    //int cutoffScore = 100;
                    //double eValueThreshold = 10.0;
                    
                    //if (bitScore >= cutoffScore && eValue < eValueThreshold) {
                        //fprintf(stdout, "%s,%s,%s, %f,%f,%f,%f, %d,%d,%d,%f,%d\n", 
                        //fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d,%d\n", 
                            //qid.c_str(), qHeader.substr(1,1000).c_str(), seqId.c_str(), 
                            ////pIdentityGapped, pIdentityUngapped, pIdentityGapOpeningOnly, pCoverage,
                            //alignLen, qStart, qEnd, 
                            //eValue, genericScore, bitScore);
                    //}
                    if (!check_exclusion(qGi, seqId, qCutLocStart, qCutLocEnd, sStart, sEnd, 100)) {
                        
                        char blastRes[MAXSTR];
                        
                        /// DEBUG
                        fprintf(stdout, "%s,%s,%s,%d,%d,%d,%1.e,%d,%d\n", 
                            qid.c_str(), qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            alignLen, qStart, qEnd, 
                            eValue, genericScore, bitScore);
                        
                        sprintf(blastRes, "%s,%s,%d,%d,%d,%1.e,%d,%d", 
                            qHeader.substr(1,MAXSTR).c_str(), seqId.c_str(), 
                            alignLen, qStart, qEnd, 
                            eValue, genericScore, bitScore);
                        
                        char *newKey = (char*)(qid.c_str());

                        /// ADD <KEY = "QUERYID,CHUNKNAME", VALUE="BLASTvecResULT"> TO KV
                        kv->add(newKey, strlen(newKey) + 1, blastRes, strlen(blastRes) + 1);
                    }
                    else 
                        fprintf(stderr, "### INFO: EXCLUSION - qGi=%s sGi=%s \
                            qCutLocStart=%d, qCutLocEnd=%d, sStart=%d, sEnd=%d\n",
                            qGi.c_str(), seqId.c_str(), qCutLocStart, qCutLocEnd, sStart, sEnd);
                        
                }
            }
            //*/
        }
    } // END WHILE
}

/* ---------------------------------------------------------------------- */
bool check_exclusion(string qGi, string sGi, int qCutLocStart, 
    size_t qCutLocEnd, size_t sStart, size_t sEnd, 
    size_t threshold)
/* ---------------------------------------------------------------------- */
{   
    ///
    /// To exclude Blast result from the original sequence from which 
    /// the input query is originated (sampled). Basically if qGi == sGi 
    /// and qCutLocStart is similar with sStart and qCutLocEnd is similar 
    /// with sEnd in terms of coordinates, the result should be excluded.
    /// 
    /// Orig seq: ----------------XXXXXXXXXXXXXXXX----------------------
    ///                           |              |
    ///                    qCutLocStart      qCutLocEnd 
    ///
    /// Query:                    XXXXXXXXXXXXXXXX
    ///                             |          | 
    ///                          qStart       qEnd
    ///
    /// Subject:  ------------------XXXXXXXXXXXX------------------------
    ///                             |          |
    ///                           sStart      sEnd
    ///
    bool ret = false;
    
    if (qGi == sGi) {
        if (qCutLocStart < 0) {            
            /// In >gi|222299657|18|3605|-400|3604
            /// -400|3604 means query[-400:3604] in Python.
            qCutLocStart = qCutLocEnd + 1 - qCutLocStart;
            qCutLocEnd += 1;
        }        
        if ((qCutLocStart - threshold <= sStart && sStart <= qCutLocStart + threshold) &&
            (qCutLocEnd - threshold <= sEnd && sEnd <= qCutLocEnd + threshold)) 
            ret = true;
    }    
    
    return ret;
}

/* ---------------------------------------------------------------------- */
void set_default_opts(CRef<CBlastOptionsHandle> optsHandle)
/* ---------------------------------------------------------------------- */
{
    optsHandle->SetEvalueThreshold(10);
    //optsHandle->SetMatchReward(0);
    //optsHandle->SetMismatchPenalty(0);
    //optsHandle->SetMatrixName("BLOSUM62");
    optsHandle->SetCutoffScore(100);
    if (CBlastNucleotideOptionsHandle* nuclHandle =
            dynamic_cast<CBlastNucleotideOptionsHandle*>(&*optsHandle)) {

        //nuclHandle->SetMatchReward(0);
        //nuclHandle->SetMismatchPenalty(0);
        nuclHandle->SetMatrixName("BLOSUM62");
        nuclHandle->SetCutoffScore(100);
        nuclHandle->SetEvalueThreshold(10);
    }

    return;
}

///* ---------------------------------------------------------------------- */
//int key_compare_str(char *p1, int len1, char *p2, int len2)
///* ---------------------------------------------------------------------- */
//{
    ////int i1 = *(int *) p1;
    ////int i2 = *(int *) p2;
    //return strcmp(p1, p2);
    ////if (i1 > i2) return -1;
    ////else if (i1 < i2) return 1;
    ////else return 0;
//}

/* ---------------------------------------------------------------------- */
int key_compare_int(char *p1, int len1, char *p2, int len2)
/* ---------------------------------------------------------------------- */
{
    //int i1 = *(int *) p1;
    //int i2 = *(int *) p2;
    size_t i1 = str2uint(string(p1));
    size_t i2 = str2uint(string(p2));
    /// Ascending order
    if (i1 > i2) return 1;
    else if (i1 < i2) return -1;
    /// Descending order
    //if (i1 > i2) return -1;
    //else if (i1 < i2) return 1;
    else return 0;
}

///* ---------------------------------------------------------------------- */
//void get_node_name(int itask, KeyValue *kv, void *ptr)
///* ---------------------------------------------------------------------- */
//{
    //procName *nn = (procName *) ptr;
    //char id[MAXSTR];
    //sprintf(id, "%d", nn->myId);
    //kv->add(nn->pName, strlen(nn->pName) + 1, id, strlen(id) + 1);
//}

/* ---------------------------------------------------------------------- */
//void collect_node_names(char *key, int keybytes, char *multivalue,
//                        int nvalues, int *valuebytes, KeyValue *kv,
//                        void *ptr)
/* ---------------------------------------------------------------------- */
//{
//kv->add(key, keybytes, multivalue, sizeof(multivalue));
//}

///* ---------------------------------------------------------------------- */
//void save_node_names(uint64_t itask, char *key, int keybytes, char *value,
                     //int valuebytes, KeyValue *kv, void *ptr)
///* ---------------------------------------------------------------------- */
//{
    ////multimap<string, unsigned int> *nNames = (multimap<string, unsigned int> *) ptr;
    //vector<string>* vecNames = (vector<string> *) ptr;
    //string pName(key);
    ////size_t myid = atoi(value);
    ////nNames->insert(pair<string, unsigned int>(pName, myid));
    //vector<string>::const_iterator loc = find(vecNames->begin(), vecNames->end(), pName);
    //if (loc == vecNames->end())
        //vecNames->push_back(pName);
//}


/// EOF

