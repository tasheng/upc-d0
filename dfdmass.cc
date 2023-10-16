#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>
#include <TMath.h>
#include <TLatex.h>
#include <TSystem.h>

using namespace ROOT;


void dfdmass() {
  ROOT::EnableImplicitMT(); // Enable ROOT's implicit multi-threading

  const char *input = "/home/data/run3RapidValidation/"
                      "run374719_374728_374729_374730_374731.root";
  TString dir = "/home/data/run3RapidValidation/";
  // TFile *finput = TFile::Open(input);
  // TTree *nt = (TTree *)finput->Get("Dfinder/ntDkpi");
  TChain *nt = new TChain("Dfinder/ntDkpi");
  nt->Add(dir + "*lowerDcut.root");

  TChain *hlt_hlt = new TChain("hltanalysis/HltTree");
  TChain *skim_hlt = new TChain("skimanalysis/HltTree");
  TChain *zdc = new TChain("zdcanalyzer/zdcdigi");
  TChain *hi = new TChain("hiEvtAnalyzer/HiTree");
  TChain *track = new TChain("ppTracks/trackTree");
  TChain *jet = new TChain("ak4CaloJetAnalyzer/caloJetTree"); // jtpt,jtphi
  std::vector<TChain *> friends = {hlt_hlt, skim_hlt, zdc, hi, track, jet};
  for (auto f : friends) {
    f->Add(dir + "*lowerDcut.root");
    nt->AddFriend(f);
  }

  RDataFrame df(*nt);
  auto entries = df.Count();
  std::cout << "entries: " << *entries << "\n";

  auto df_cut =
      df.Filter("pprimaryVertexFilter && !HLT_HIMinimumBiasHF1ANDZDC1nOR_v1"
                "&& TMath::Abs(PVz)<15")
          .Define("pmcut",
                  "(sumPlus>1100&&sumMinus<1100||sumPlus<1100&&sumMinus>1100)")
          .Define("trketacut", "abs(Dtrk1Eta)<2.4"
                               "&& abs(Dtrk2Eta)<2.4")
          .Define("trkNHitcut", "(Dtrk1PixelHit+Dtrk1StripHit)>=11"
                                "&& (Dtrk2PixelHit+Dtrk2StripHit)>=11")
          .Define("trkPtErrcut", "abs(Dtrk1PtErr/Dtrk1Pt)<0.3"
                                 "&& abs(Dtrk2PtErr/Dtrk2Pt)<0.3")
          // .Define("ycut", "abs(Dy) < 1")
          // .Define("ycut", "abs(Dy) > 1 && abs(Dy) < 2")
          .Define("ptcut", "Dpt>2")
          .Define("alphacut", "Dalpha<0.15")
          .Define("chi2cut", "Dchi2cl>0.05")
          .Define("dlscut", "DsvpvDisErr<10")
          .Define("dthetacut", "cos(Ddtheta) > 0.8")
          .Define("lxybscut", "DlxyBS / DlxyBSErr > 1")
          // .Define("allcut",
          //         "pmcut && trketacut && trkNHitcut && trkPtErrcut &&"
          //         "ycut && ptcut && alphacut && chi2cut && dlscut &&"
          //         "dthetacut && lxybscut")
          .Define("allcut", "pmcut && trketacut && trkNHitcut &&"
                            "chi2cut && dlscut &&"
                            "dthetacut")
          .Define("DmassGood", "Dmass[allcut]")
          .Filter("! Dmass[allcut].empty()");
  entries = df_cut.Count();
  // entries = df_cut.Sum("DmassGood");
  std::cout << *entries << "\n";

  // TH1F *hmass = new TH1F("hmass", "D^{0} mass; D^{0} mass/GeV; Entries/0.003
  // GeV", 100, 1.75, 2.05);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);
  gStyle->SetOptStat(0);
  TFile *f = new TFile("output.root", "recreate");

  TLatex lat;
  lat.SetNDC();

  // std::vector<float> ptlist = {2, 5, 8, 11, 14};
  std::vector<float> ptlist = {0, 1, 2, 4, 6, 8, 10, 12, 14, 16};
  std::vector<float> ylist = {0, 1, 1.5, 2, 2.4};
  for (auto i = 0; i < ptlist.size() - 1; ++i) {
    std::cout << "filtering " << i << "\n";
    TString filt = Form("dpt_%.0f_%.0f", ptlist[i], ptlist[i + 1]);
    auto df_ptsel = df_cut.Define(filt.Data(),
                                  Form("Dpt > %f && Dpt < %f", ptlist[i], ptlist[i + 1]))
      .Define("dmasspt", Form("Dmass[allcut && %s]", filt.Data()));
    auto hmeta2 = df_ptsel.Histo1D({Form("hmass_%d", i),
        "D^{0} mass; D^{0} mass/GeV; Entries/0.003 GeV", 100,
         1.75, 2.05},
        "dmasspt");
    hmeta2->SetLineColor(kBlack);
    hmeta2->SetLineWidth(2);
    f->cd();
    hmeta2->Write();

    canvas->cd();
    hmeta2->Draw("e");
    lat.DrawLatex(0.55, 0.7, "|y| < 1");
    // lat.DrawLatex(0.55, 0.7, "1 < |y| < 2");
    lat.DrawLatex(0.55, 0.6,
                  Form("%.0f < p_{T} < %.0f", ptlist[i], ptlist[i + 1]));
    gSystem->mkdir("png");
    gSystem->mkdir("pdf");
    // canvas->SaveAs("pdf/dmass_y_1_2_" + filt + ".pdf");
    // canvas->SaveAs("png/dmass_y_1_2_" + filt + ".png");
    canvas->SaveAs("pdf/dmass_y_0_1_" + filt + ".pdf");
    canvas->SaveAs("png/dmass_y_0_1_" + filt + ".png");
  }
}
