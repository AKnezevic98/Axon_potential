#include <math.h>
#include "pbPlots.h"
#include "supportLib.h"

/*
This program simulates potential on an axon.

The program plots axon potential through time as well as concentration of open channels for Potassium, Sodium and leaky channels.

The program also plots number of triggers of axon as a function of external potential.
*/

unsigned char main(){
	
	//relevant functions
	double an(double x){
		return 0.01*(x+55)/(1-exp(-(x+55)/10));
	}
	double bn( double x){
		return 0.125*exp(-(x+65)/80);
	}
	double am(double x){
		return 0.1*(x+40)/(1-exp(-(x+40)/10));
	}
	double bm( double x){
		return 4*exp(-(x+65)/18);
	}
	double ah(double x){
		return 0.07*exp(-(x+65)/20);
	}
	double bh( double x){
		return 1/(1+exp(-(x+35)/10));
	}
	
	double noo(double x){
		return (0.01*(x+55)/(1-exp(-(x+55)/10)))/(0.01*(x+55)/(1-exp(-(x+55)/10)) + 0.125*exp(-(x+65)/80));
	}
	double moo(double x){
		return (0.1*(x+40)/(1-exp(-(x+40)/10)))/(0.1*(x+40)/(1-exp(-(x+40)/10)) + 4*exp(-(x+65)/18));
	}
	double hoo(double x){
		return (0.07*exp(-(x+65)/20))/(0.07*exp(-(x+65)/20) + 1/(1+exp(-(x+35)/10)));
	}
	
	//more relevant functions
	double min(double a[], unsigned long int length){
		double minimum = a[0];
		unsigned long int i=0;
		for(i=0; i<length; i++){
			if(a[i]<minimum){
				minimum = a[i];
			}
		}
		return minimum;
	}
	double max(double a[], unsigned long int length){
		double maximum = a[0];
		unsigned long int i=0;
		for(i=0; i<length; i++){
			if(a[i]>maximum){
				maximum = a[i];
			}
		}
		return maximum;
	}
	
	double trigger_number(double a[], double t[], unsigned long int length){
		double num = 0;
		double last = 0;
		unsigned long int i=0;
		for(i=1; i<length-1; i++){
			if(a[i]>-55 && a[i-1]<a[i] && a[i+1]<a[i]){
				num++;
				if(t[i]-last<2.8){ //t[i+1]>=tc && t[i]<=tc && a[i]<-15 || 
					num--;
				}
				last = t[i];
			}
		}
		return num;
	}
	
	
	unsigned long i = 0;
	
	double gNa = 120, gK = 36, gL = 0.3, ENa = 50, EK = -77, EL = -54.4; //constants of the model
	
	double V0 = -65, n0 = 0.32, m0 = 0.053, h0 = 0.6; //initial conditions of the model
	
	double Vext0 = 54.6, tstart = 20, tend = 80; //parameters of the external influence
	
	double dt = 0.003; //timestep parameter
	double tmax = 100; //max t
	unsigned long N = (unsigned long)floor(tmax/dt)+1; //number of elements in arrays
	
	double ts[N], Vs[N], Vext[N], ns[N], ms[N], hs[N]; //initialising arrays
	Vs[0] = V0; ns[0] = n0; ms[0] = m0; hs[0] = h0;
	
	// filling initial arrays
	for(i=0; i<N; i++){
		ts[i] = i*dt;
		if(ts[i]>=tstart && ts[i]<=tend){
			Vext[i] = Vext0;
		} else {
			Vext[i] = 0;
		}
	}
	
	//integration of the equations
	for(i=1; i<N; i++){
		Vs[i] = Vs[i-1] - dt*( gNa*pow(ms[i-1],3)*hs[i-1]*(Vs[i-1]-ENa) + gK*pow(ns[i-1],4)*(Vs[i-1]-EK) + gL*(Vs[i-1]-EL) - Vext[i] );
		ns[i] = ns[i-1] + dt*(an(Vs[i-1]) + bn(Vs[i-1]))*(noo(Vs[i-1])-ns[i-1]);
		ms[i] = ms[i-1] + dt*(am(Vs[i-1]) + bm(Vs[i-1]))*(moo(Vs[i-1])-ms[i-1]);
		hs[i] = hs[i-1] + dt*(ah(Vs[i-1]) + bh(Vs[i-1]))*(hoo(Vs[i-1])-hs[i-1]);
	}
	
	
	//plotting part
	
	double xmin, xmax, ymin, ymax, ymin1, ymin2, ymin3, ymax1, ymax2, ymax3, cmin, cmax;
	
	xmin = min(ts, N); xmax = max(ts, N);
	ymin = min(Vs, N); ymax = max(Vs, N);
	
	ymin1 = min(ns, N); ymax1 = max(ns, N);
	ymin2 = min(ms, N); ymax2 = max(ms, N);
	ymin3 = min(hs, N); ymax3 = max(hs, N);
	
	if(ymin1<ymin2 && ymin1<ymin3){
		cmin = ymin1;
	}else if(ymin2<ymin3){
		cmin = ymin2;
	}else{
		cmin = ymin3;
	}
	if(ymax1>ymax2 && ymax1>ymax3){
		cmax = ymax1;
	}else if(ymax2>ymax3){
		cmax = ymax2;
	}else{
		cmax = ymax3;
	}
	
	
	//plotting of axion potential
	
	ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = ts;
	series->xsLength = sizeof(ts)/sizeof(double);
	series->ys = Vs;
	series->ysLength = sizeof(Vs)/sizeof(double);
	series->linearInterpolation = true;
	series->lineThickness = 2;
	series->color = CreateRGBColor(1, 0, 0);

	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 1920;
	settings->height = 1080;
	settings->autoBoundaries = false;
	settings->xMin = xmin;
	settings->xMax = xmax;
	settings->yMin = ymin;
	settings->yMax = ymax;
	settings->autoPadding = true;
	settings->title = L"Axion potential vs time";
	settings->titleLength = wcslen(settings->title);
	settings->xLabelLength = wcslen(settings->xLabel);
	settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s [] = {series};
	settings->scatterPlotSeries = s;
	settings->scatterPlotSeriesLength = 1;

	RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	size_t length;
	double *pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, "Axion_potential.png");
	
	
	//plotting of n, m and h
	
	series = GetDefaultScatterPlotSeriesSettings();
	series->xs = ts;
	series->xsLength = sizeof(ts)/sizeof(double);
	series->ys = ns;
	series->ysLength = sizeof(ns)/sizeof(double);
	series->linearInterpolation = true;
	series->lineThickness = 2;
	series->color = CreateRGBColor(1, 0, 0);
	
	ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = ts;
	series2->xsLength = sizeof(ts)/sizeof(double);
	series2->ys = ms;
	series2->ysLength = sizeof(ms)/sizeof(double);
	series2->linearInterpolation = true;
	series2->lineThickness = 2;
	series2->color = CreateRGBColor(0, 1, 0);
	
	ScatterPlotSeries *series3 = GetDefaultScatterPlotSeriesSettings();
	series3->xs = ts;
	series3->xsLength = sizeof(ts)/sizeof(double);
	series3->ys = hs;
	series3->ysLength = sizeof(hs)/sizeof(double);
	series3->linearInterpolation = true;
	series3->lineThickness = 2;
	series3->color = CreateRGBColor(0, 0, 1);
	
	settings = GetDefaultScatterPlotSettings();
	settings->width = 1920;
	settings->height = 1080;
	settings->autoBoundaries = false;
	settings->xMin = xmin;
	settings->xMax = xmax;
	settings->yMin = cmin;
	settings->yMax = cmax;
	settings->autoPadding = true;
	settings->title = L"n(red), m(green) and h(blue) vs time";
	settings->titleLength = wcslen(settings->title);
	settings->xLabelLength = wcslen(settings->xLabel);
	settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s1 [] = {series, series2, series3};
	settings->scatterPlotSeries = s1;
	settings->scatterPlotSeriesLength = 3;

	canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, "NMH.png");
	
	free(series2); free(series3);
	
	
	
	//plotting number of number of triggers as a function of strength of external influence
	
	double Vmax = 185, dV = 0.1;
	unsigned long n = (unsigned long)floor(Vmax/dV)+1;
	unsigned long j = 0;
	double nums[n];
	double Vexts[n];
	
	for(j=0; j<n; j++){
		Vexts[j] = j*dV;
		Vext0 = Vexts[j];
		
		for(i=0; i<N; i++){
			if(ts[i]>=tstart && ts[i]<=tend){
				Vext[i] = Vext0;
			} else {
				Vext[i] = 0;
			}
		}
		
		//integration of the equations
		for(i=1; i<N; i++){
			Vs[i] = Vs[i-1] - dt*( gNa*pow(ms[i-1],3)*hs[i-1]*(Vs[i-1]-ENa) + gK*pow(ns[i-1],4)*(Vs[i-1]-EK) + gL*(Vs[i-1]-EL) ) + dt*Vext[i];
			ns[i] = ns[i-1] + dt*(an(Vs[i-1]) + bn(Vs[i-1]))*(noo(Vs[i-1])-ns[i-1]);
			ms[i] = ms[i-1] + dt*(am(Vs[i-1]) + bm(Vs[i-1]))*(moo(Vs[i-1])-ms[i-1]);
			hs[i] = hs[i-1] + dt*(ah(Vs[i-1]) + bh(Vs[i-1]))*(hoo(Vs[i-1])-hs[i-1]);
		}
		
		nums[j] = trigger_number(Vs, ts, N);
	}
	
	
	xmin = min(Vexts, n); xmax = max(Vexts, n);
	ymin = min(nums, n)-1; ymax = max(nums, n)+1;
	
	series = GetDefaultScatterPlotSeriesSettings();
	series->xs = Vexts;
	series->xsLength = sizeof(Vexts)/sizeof(double);
	series->ys = nums;
	series->ysLength = sizeof(nums)/sizeof(double);
	series->linearInterpolation = true;
	series->lineThickness = 2;
	series->color = CreateRGBColor(0, 0, 1);

	settings = GetDefaultScatterPlotSettings();
	settings->width = 1920;
	settings->height = 1080;
	settings->autoBoundaries = false;
	settings->xMin = xmin;
	settings->xMax = xmax;
	settings->yMin = ymin;
	settings->yMax = ymax;
	settings->autoPadding = true;
	settings->title = L"Number of triggers as a function of Vext";
	settings->titleLength = wcslen(settings->title);
	settings->xLabelLength = wcslen(settings->xLabel);
	settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s2 [] = {series};
	settings->scatterPlotSeries = s2;
	settings->scatterPlotSeriesLength = 1;

	canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, "number_of_triggers.png");
	
	DeleteImage(canvasReference->image);
	free(series);
	free(settings);
	free(pngdata);
	free(canvasReference);
	
	
	return 0;
}
