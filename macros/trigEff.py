#!/usr/bin/env python

import os, sys
import ROOT

def getEff (fIn,hallIn,hpassIn): 
  
  print fIn, " ", hallIn, " ", hpassIn

  fIn = ROOT.TFile(fIn,"READ")
  
  hall = fIn.Get("trig/"+hallIn)
  hpass = fIn.Get("trig/"+hpassIn)
  
  greff = ROOT.TGraphAsymmErrors(hpass, hall, "cp")
  greff.SetName("eff_"+hallIn+hpassIn.lstrip(hallIn)) 

  fout = ROOT.TFile(greff.GetName()+".root","RECREATE") 
  fout.cd()
  greff.Write()
  fout.Close() 
  
  ROOT.gROOT.SetStyle('Plain')
  ROOT.gStyle.SetOptStat(0)

  ceff = ROOT.TCanvas("eff_"+hallIn+hpassIn.lstrip(hallIn),"eff_"+hallIn+hpassIn.lstrip(hallIn),800,600)
  ceff.cd()
  greff.Draw("ap")
  ceff.SaveAs(ceff.GetName()+".pdf")
  ceff.SaveAs(ceff.GetName()+".C")
  
  fIn.Close()

getEff("triggers_data_Run2015D_05Oct2015.root","htak4","htak4_HLT_PFHT800_v")
getEff("triggers_data_Run2015D_05Oct2015.root","htak4","htak4_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")
getEff("triggers_data_Run2015D_05Oct2015.root","htak4","htak4_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_v")

getEff("triggers_data_Run2015D_05Oct2015.root","htak8","htak8_HLT_PFHT800_v")
getEff("triggers_data_Run2015D_05Oct2015.root","htak8","htak8_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")
getEff("triggers_data_Run2015D_05Oct2015.root","htak8","htak8_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_v")

getEff("triggers_data_Run2015D_05Oct2015.root","ptak82nd","ptak82nd_HLT_PFHT800_v")
getEff("triggers_data_Run2015D_05Oct2015.root","ptak82nd","ptak82nd_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")
getEff("triggers_data_Run2015D_05Oct2015.root","ptak82nd","ptak82nd_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_v")
