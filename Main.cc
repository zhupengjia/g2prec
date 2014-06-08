// -*- C++ -*-

/* Main.cc
 * Main function of this simulation.
 * It parses command line parameters use getopt (GNU C Library).
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Mar 2013, C. Gu, Use Run.C script to do simulation
//   Sep 2013, C. Gu, Since G2PRun will take care of the configuration files, make it simpler here.
//   Feb 2014, C. Gu, Modified for reconstruction.
//

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <ctype.h>
#include <getopt.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "THaEvent.h"

#include "libconfig.h++"

#include "G2PRec.hh"

//#define DEBUG

using namespace std;
using namespace libconfig;

Config* gConfig = new Config();
G2PRec* gRec = NULL;

const char* DBDIR = "./recdb";
const char* ROOTDIR = ".";

bool isexist(const char* filename);
void usage(int argc, char** argv);

int Configure(int run)
{
    static const char* const here = "Main::Configure()";

    TFile *f = new TFile(Form("%s/g2p_%d.root", ROOTDIR, run), "READ");
    TTree *t = (TTree *) f->Get("T");

    THaEvent *event = new THaEvent();
    t->SetBranchAddress("Event_Branch", &event);

    t->GetEntry(0);

    int label = int(event->GetHeader()->GetEvtTime() / 1.0e6);

    delete event;
    event = NULL;

    f->Close();

    void* dirp = gSystem->OpenDirectory(DBDIR);
    if (dirp == NULL) return -1;

    const char* result;
    string item;
    vector<int> time;
    while ((result = gSystem->GetDirEntry(dirp)) != NULL) {
        item = result;
        size_t pos;
        for (pos = 0; pos < item.length(); pos++) {
            if (!isdigit(item[pos])) break;
        }
        if (pos == item.length()) {
            time.push_back(atoi(item.c_str()));
        }
    }
    gSystem->FreeDirectory(dirp);

    if (time.size() > 0) {
        sort(time.begin(), time.end());
        vector<int>::iterator it;
        for (it = time.begin(); it != time.end(); it++) {
            if (label < (*it)) {
                if (it == time.begin()) return -1;
                it--;
                break;
            }
        }
        if (it == time.end()) it--; // Assume the last directory is always valid

        // Only work for g2p
        if (run < 20000) {
            Info(here, "Use %s", Form("%s/%d/db_L.optics.cfg", DBDIR, (*it)));
            gConfig->readFile(Form("%s/%d/db_L.optics.cfg", DBDIR, (*it)));
        } else if (run < 40000) {
            Info(here, "Use %s", Form("%s/%d/db_R.optics.cfg", DBDIR, (*it)));
            gConfig->readFile(Form("%s/%d/db_R.optics.cfg", DBDIR, (*it)));
        } else return -1;
    }

    gRec = new G2PRec();

    const char* rundbfile = Form("%s/db_rec.dat", DBDIR);
    if (!isexist(rundbfile)) {
        return -1;
    }

    FILE *fp = fopen(rundbfile, "r");

    int id = 0;
    double e = 0.0, p = 0.0;
    bool found = false;
    while (!feof(fp)) {
        fscanf(fp, "%d%le%le", &id, &e, &p);
        if (id == run) {
            gRec->SetBeamEnergy(e);
            gRec->SetHRSMomentum(p);
            found = true;
            break;
        }
    }
    if (!found) return -1;

    fclose(fp);

    return 0;
}

int Insert(int run)
{
    static const char* const here = "Main::Insert()";

    int fDebug;
    gConfig->lookupValue("debug", fDebug);

    const char* arm = "R";
    if (run < 20000) arm = "L";

    const char* filename = Form("%s/g2p_%d.root", ROOTDIR, run);
    int inc = 0;

    while (isexist(filename)) {
        TFile *f = new TFile(filename, "UPDATE");
        Info(here, "Opening existed rootfile %s ...", filename);

        TTree *t = (TTree *) f->Get("T");

        float fV5bpm_bpm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        float fV5bpmave_bpm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        int fBPMAvail;
        double fV5tp_tr[5];
        THaEvent *event = new THaEvent();

        t->SetBranchAddress("Event_Branch", &event);

        if (t->FindBranch(Form("%srb.tgt_m13_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m13_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m13_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m13_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m13_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m13_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m13_y", arm), &fV5bpmave_bpm[2]);
            fV5bpm_bpm[4] = -13.6271;
        } else if (t->FindBranch(Form("%srb.tgt_m12_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m12_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m12_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m12_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m12_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m12_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m12_y", arm), &fV5bpmave_bpm[2]);
            fV5bpm_bpm[4] = -12.5476;
        } else if (t->FindBranch(Form("%srb.tgt_m10_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m10_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m10_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m10_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m10_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m10_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m10_y", arm), &fV5bpmave_bpm[2]);
            fV5bpm_bpm[4] = -10.81;
        } else {
            t->SetBranchAddress(Form("%srb.tgt_0_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_0_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_0_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_0_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_0_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_0_y", arm), &fV5bpmave_bpm[2]);
            fV5bpm_bpm[4] = 0.0;
        }

        t->SetBranchAddress(Form("%srb.bpmavail", arm), &fBPMAvail);

        t->SetBranchAddress(Form("%s.gold.x", arm), &fV5tp_tr[0]);
        t->SetBranchAddress(Form("%s.gold.th", arm), &fV5tp_tr[1]);
        t->SetBranchAddress(Form("%s.gold.y", arm), &fV5tp_tr[2]);
        t->SetBranchAddress(Form("%s.gold.ph", arm), &fV5tp_tr[3]);
        t->SetBranchAddress(Form("%s.gold.dp", arm), &fV5tp_tr[4]);

        double fV5rec_tr[5];
        double fV5rec_lab[5];
        double fV5corr_tr[5];

        TList newBranch;

        newBranch.Add(t->Branch(Form("%s.rec.x", arm), &fV5rec_tr[0], Form("%s.rec.x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.th", arm), &fV5rec_tr[1], Form("%s.rec.th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.y", arm), &fV5rec_tr[2], Form("%s.rec.y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.ph", arm), &fV5rec_tr[3], Form("%s.rec.ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.dp", arm), &fV5rec_tr[4], Form("%s.rec.dp/D", arm)));

        newBranch.Add(t->Branch(Form("%s.rec.l_x", arm), &fV5rec_lab[0], Form("%s.rec.l_x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_th", arm), &fV5rec_lab[1], Form("%s.rec.l_th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_y", arm), &fV5rec_lab[2], Form("%s.rec.l_y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_ph", arm), &fV5rec_lab[3], Form("%s.rec.l_ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_z", arm), &fV5rec_lab[4], Form("%s.rec.l_z/D", arm)));

#ifdef DEBUG
        newBranch.Add(t->Branch(Form("%s.corr.x", arm), &fV5corr_tr[0], Form("%s.corr.x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.th", arm), &fV5corr_tr[1], Form("%s.corr.th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.y", arm), &fV5corr_tr[2], Form("%s.corr.y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.ph", arm), &fV5corr_tr[3], Form("%s.corr.ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.dp", arm), &fV5corr_tr[4], Form("%s.corr.dp/D", arm)));
#endif

        int N = t->GetEntries();
        int evnum;
        TIter next(&newBranch);
        for (int i = 0; i < N; i++) {
            t->GetEntry(i);
            evnum = Int_t(event->GetHeader()->GetEvtNum());
            if (fDebug > 0) {
                Info(here, "Processing event %d ......", evnum);
            } else if ((i % 10000 == 0)&&(i != 0)) {
                Info(here, "%d events Processed ......", i);
            }
            if ((fBPMAvail < 0.5) || (fV5tp_tr[0] > 1.0e8) || (fV5tp_tr[1] > 1.0e8) || (fV5tp_tr[2] > 1.0e8)
                    || (fV5tp_tr[3] > 1.0e8) || (fV5tp_tr[4] > 1.0e8)) {
                for (int i = 0; i < 5; i++) {
                    fV5rec_tr[i] = 1e+38;
                    fV5rec_lab[i] = 1e+38;
                }
            } else {
                gRec->Process(fV5bpm_bpm, fV5bpmave_bpm, fV5tp_tr, fV5corr_tr, fV5rec_tr, fV5rec_lab);
            }
            next.Reset();
            while (TBranch * br = (TBranch*) next()) {
                br->Fill();
            }
        }

        t->Write("", TObject::kOverwrite);
        f->Close();

        inc++;
        filename = Form("%s/g2p_%d_%d.root", ROOTDIR, run, inc);
    }

    Info(here, "No more rootfiles.");

    return 0;
}

int main(int argc, char** argv)
{
    int c;

    int fRun = 0;

    while (1) {
        static struct option long_options[] = {
            {"dbdir", required_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"rootdir", required_argument, 0, 'r'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "d:hr:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 0:
            break;
        case 'd':
            DBDIR = Form("%s", optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
        case 'r':
            ROOTDIR = Form("%s", optarg);
            break;
        case '?':
        default:
            usage(argc, argv);
            exit(-1);
        }
    }

    if (optind < argc) {
        fRun = atoi(argv[optind++]);
    } else {
        usage(argc, argv);
        exit(-1);
    }

    const char* fFileName = Form("%s/g2p_%d.root", ROOTDIR, fRun);

    if (!isexist(fFileName)) {
        usage(argc, argv);
        exit(-1);
    }

    if (Configure(fRun) != 0) {
        Error("Main()", "Please check database!");
        exit(-1);
    }

    Insert(fRun);

    delete gRec;
    gRec = NULL;

    return 0;
}

bool isexist(const char* fname)
{
    return (access(fname, F_OK) == 0) ? true : false;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] RunNumber \n", argv[0]);
    printf("  -d, --dbdir=$PWD/recdb         Set db directory\n");
    printf("  -h, --help                     Print this small usage guide\n");
    printf("  -r, --rootdir=$PWD             Set db directory\n");
}