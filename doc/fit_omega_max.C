#include <cmath>

double F_E1(double w) {
   double result;

   result = tanh(w/2.);
   return(0.25 * (1. - result));
}

double F_E2(double w) {
   double result;

   result = sinh(w)*(cosh(w) + 2.) / pow(cosh(w) + 1., 2);
   return(0.25 * (1. - result));
}

double F_E3(double w) {
   double result;

   result = sinh(w)*(6*cosh(w)+cosh(2*w)+8.);
   result /= 2 * pow(cosh(w) + 1., 3);
   return(0.25 * (1. - result));
}

double F_E4(double w) {
   double result;

   result = 56 * sinh(w);
   result += 28 * sinh(2*w);
   result += 8 * sinh(3*w);
   result += sinh(4*w);
   result /= 8 * pow(cosh(w) + 1, 4);
   return(0.25 * (1. - result));
}

double F_E5(double w) {
   double result;

   result = 130 * cosh(w);
   result += 46 * cosh(2*w);
   result += 10 * cosh(3*w);
   result += cosh(4*w);
   result += 128;
   result *= sinh(w);
   result /= 8 * pow(cosh(w) + 1, 5);
   return(0.25 * (1. - result));
}

double F_E6(double w) {
   double result;

   result = 792 * sinh(w);
   result += 495 * sinh(2*w);
   result += 220 * sinh(3*w);
   result += 66 * sinh(4*w);
   result += 12 * sinh(5*w);
   result += sinh(6*w);
   result /= 32 * pow(cosh(w) + 1, 6);
   return(0.25 * (1. - result));
}

double real_func(int lambda, double w) {
   switch(lambda) {
    case 1:
      return(F_E1(w));
      break;
    case 2:
      return(F_E2(w));
      break;
    case 3:
      return(F_E3(w));
      break;
    case 4:
      return(F_E4(w));
      break;
    case 5:
      return(F_E5(w));
      break;
    case 6:
      return(F_E6(w));
      break;
   }
   return(-1);
}

double fit_func(int lambda, double a, double w) {

   double result;

   result = exp((double)lambda * (a - w));
   return(result);
}

void fit_omega_max() {

   double w;
   double wm[] = {10.82, 5.96, 4.37, 3.59, 3.13, 2.88};
   double a[] = {-.693, .203, .536, .716, .829, .962};
   double a2[] = {-.69317, .19915, .51685, .67030, .75467, .79585};
   char line[1024];
   double dw = 0.0005, first, y;
   
   TCanvas *c1 = new TCanvas("c1", "Gosia omega limits", 1024, 768);
   c1->Divide(3,2);

   // Loop on multipolarities
   for (int lambda = 1; lambda <= 6; lambda++) {

      // Create a graph
      TGraph *g = new TGraph(1001);
      TGraph *g2 = new TGraph(1001);
      sprintf(line, "E%d", lambda);
      g->SetTitle(line);
      first = 0;
      for (int i = 0; i <= 1000; i++) {
	 w = wm[lambda - 1] + (double)(i - 500) * dw;
	 y = log(real_func(lambda, w));
	 if (first == 0) first = w;
	 g->SetPoint(i, y, w);
	 g2->SetPoint(i, log(fit_func(lambda, a[lambda-1], w)), w);

      }

      // Perform a fit
      w = a2[lambda-1] - log(1e-5) / (double)lambda;
      TF1 *fit = new TF1("fit", "pol1(0)", w - 500 * dw, w + 500 * dw);
      fit->SetLineColor(2);
      fit->SetLineWidth(0.008);
      fit->SetParameter(1, -1./(double)lambda);
      g->Fit(fit);
      w= (fit->GetParameter(0) - log(1e-5)) / (double)lambda;
      sprintf(line,"%.5lf -1/%.5lf * ln(a_c)\n", fit->GetParameter(0), -1. / fit->GetParameter(1));
      delete fit;

      // Select pad
      c1->cd(lambda);

      // Setup graph
      g->GetXaxis()->SetTitle("ln(a_c)");
      g->GetYaxis()->SetTitle("omega");
      g->GetYaxis()->SetRangeUser(wm[lambda-1]-500*dw, wm[lambda-1]+500*dw);

      // Draw graph
      g->Draw("AL");
      g2->SetLineColor(3);
      g2->Draw("same");

      // Add text to label
      TText *t1 = new TText(y, first + 0.08, "Gosia");
      TText *t2 = new TText(y, first + 0.05, "Analytic");
      TText *t3 = new TText(y, first + 0.02, line);
      t1->SetTextColor(3);
      t2->SetTextColor(1);
      t3->SetTextColor(2);
      t1->Draw("same");
      t2->Draw("same");
      t3->Draw("same");
   }
   c1->SaveAs("integration_range.eps");
   
}
