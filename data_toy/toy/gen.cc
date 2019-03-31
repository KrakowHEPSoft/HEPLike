void gen()
{
  TRandom3 *rand= new TRandom3(123);
  TFile *f = new TFile("data.root", "RECREATE");
  TTree *t = new TTree("t", "t");
  const double mean=1;
  const double sigma=2.;

  double x,w;
  t->Branch("x", &x, "x/D");
  t->Branch("w", &w, "w/D"); 

  
  for(int i=0; i< 1000; ++i)
    {
      x=rand->Gaus(mean, sigma);
      w=1.;
      t->Fill();

    }
  t->Write();




}
  
