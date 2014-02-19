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

using namespace std;
using namespace libconfig;

Config* gConfig = new Config();

bool isexist(const char* filename);
void usage(int argc, char** argv);

int Configure(int run, const char* dbdir, int debug)
{
    static const char* const here = "Main::Configure()";

    TFile *f = new TFile(Form("g2p_%d.root", run), "READ");
    TTree *t = (TTree *) f->Get("T");

    THaEvent *event = new THaEvent();
    t->SetBranchAddress("Event_Branch", &event);

    t->GetEntry(0);

    int label = int(event->GetHeader()->GetEvtTime() / 1.0e6);

    void* dirp = gSystem->OpenDirectory(dbdir);
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
            Info(here, "Use %s", Form("%s/%d/db_L.optics.cfg", dbdir, (*it)));
            gConfig->readFile(Form("%s/%d/db_L.optics.cfg", dbdir, (*it)));
        } else if (run < 40000) {
            Info(here, "Use %s", Form("%s/%d/db_R.optics.cfg", dbdir, (*it)));
            gConfig->readFile(Form("%s/%d/db_R.optics.cfg", dbdir, (*it)));
        } else return -1;
    }



    return 0;
}

int Insert()
{
    G2PRec * run = new G2PRec();

    return 0;
}

int main(int argc, char** argv)
{
    int c;

    int fDebug = 0;
    const char* fDBDir = "./recdb";
    int fRun = 0;

    while (1) {
        static struct option long_options[] = {
            {"dbdir", required_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"level", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "d:hl:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 0:
            break;
        case 'd':
            fDBDir = Form("%s", optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
        case 'l':
            fDebug = atoi(optarg);
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

    const char* fFileName = Form("g2p_%d.root", fRun);

    if (!isexist(fFileName)) {
        usage(argc, argv);
        exit(-1);
    }

    if (Configure(fRun, fDBDir, fDebug) != 0) {
        printf("Please check the database!\n");
        exit(-1);
    }

    Insert();

    return 0;
}

bool isexist(const char* fname)
{
    return (access(fname, F_OK) == 0) ? true : false;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] RunNumber \n", argv[0]);
    printf("  -d, --dbdir=./recdb            Set db directory\n");
    printf("  -h, --help                     Print this small usage guide\n");
    printf("  -l, --level=0                  Set debug level\n");
}