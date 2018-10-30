TGraphPolar * hist2tgraphPolar(TH1 * hist){

    int N = hist->GetNbinsX();

    std::vector<double> theta(N,0);
    std::vector<double> radius(N,0);
    std::vector<double> etheta(N,0);
    std::vector<double> eradius(N,0);

    for (int i=0; i<N; i++) {
        theta[i]   = (i+1)*(TMath::Pi()/(0.5*N));
        radius[i]  = hist->GetBinContent(i+1);
        etheta[i]  = TMath::Pi()/(double)N;
        eradius[i] = 0.0;
    }

    TGraphPolar * ans = new TGraphPolar(N, &theta[0], &radius[0], &etheta[0], &eradius[0]);
    return ans;

}

double mygaus(double x, double A, double mean, double sigma){

    return A*exp(-(x-mean)*(x-mean)/(2.0*sigma*sigma));
}

double fitf(Double_t *x,Double_t *par) {
        
        return  mygaus(*x, par[0],par[1],par[2])+ mygaus(*x, par[3],par[4],par[5]);
    
    }



void polar(){

    TRandom3 *rangen  = new TRandom3(0);
    int N=20;

    TH1D * test = new TH1D("polar_test","polar_test",N,-3.14159,3.1419);
    for(int i=0; i<N; i++){
        test->SetBinContent(i+1,rangen->Uniform(2,20));
    }

    TGraphPolar * grP1 = hist2tgraphPolar(test);  



    TCanvas * CPol = new TCanvas("CPol","TGraphPolar Example",500,500);
    grP1->SetTitle("");
    grP1->SetMarkerStyle(20);
    grP1->SetMarkerSize(2.);
    grP1->SetMarkerColor(4);
    grP1->SetLineColor(2);
    grP1->SetLineWidth(3);
    grP1->Draw("PE");

    TF1 *func = new TF1("fit",fitf,-3.14,3.14,6);
    func->SetParNames("Norm1","Mean1","sigma1","Norm2","Mean2","Sigma2");


    TFitResultPtr r = grP1->Fit("fit","S");

    double chi2   = r->Chi2();                  // to retrieve the fit chi2
    double par0   = r->Parameter(0);            // retrieve the value for the parameter 0
    double err0   = r->ParError(0);             // retrieve the error for the parameter 0

    // Update, otherwise GetPolargram returns 0
    CPol->Update();
    grP1->GetPolargram()->SetToRadian();

    return 0;




}
