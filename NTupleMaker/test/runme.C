void runme(){
    gROOT->ProcessLine(".L analyzer.C");
    gROOT->ProcessLine("analyzer t");
    gROOT->ProcessLine("t.Loop()");
}



