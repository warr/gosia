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
   char line[1024];
   
   TCanvas *c1 = new TCanvas("c1", "Gosia omega limits", 1024, 768);
   c1->Divide(3,2);
   
   for (int lambda = 1; lambda <= 6; lambda++) {
      TGraph *g = new TGraph(1001);
      TGraph *g2 = new TGraph(1001);
      for (int i = 0; i <= 1000; i++) {
	 w = wm[lambda - 1] + (double)(i - 500) * 0.001;
	 g->SetPoint(i, w, real_func(lambda, w));
	 g2->SetPoint(i, w, fit_func(lambda, a[lambda-1], w));
      }

      c1->cd(lambda);
      c1->GetPad(lambda)->SetLogy(true);
      g->GetYaxis()->SetRangeUser(5e-6,2e-5);
      g2->GetYaxis()->SetRangeUser(5e-6,2e-5);
      g->GetXaxis()->SetTitle("omega");
      g->GetYaxis()->SetTitle("a_c");
      g->GetXaxis()->SetRangeUser(wm[lambda-1]-.5, wm[lambda-1]+.5);
      sprintf(line, "E%d", lambda);
      g->SetTitle(line);
      g->Draw("AL");
      g2->SetLineColor(2);
      g2->Draw("same");
      TLine *l = new TLine(wm[lambda-1]-.5, 1e-5, wm[lambda-1]+.5,1e-5);
      l->Draw("same");
      TText *t1 = new TText(wm[lambda-1]+.1, 1.5e-5, "Gosia");
      TText *t2 = new TText(wm[lambda-1]+.1, 1.4e-5, "Analytic");
      t2->SetTextColor(2);
      t1->Draw("same");
      t2->Draw("same");
   }
   c1->SaveAs("integration_range.eps");
   
}
