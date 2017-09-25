#include "Main.h"
#include "MainFrame.h"
#include <ApplicationServices/ApplicationServices.h>

// called OnInit in "Main" class
IMPLEMENT_APP(Main)

bool Main::OnInit(){
	// MainFrame: public wxFrame

	ProcessSerialNumber PSN;
	GetCurrentProcess(&PSN);
	TransformProcessType(&PSN,kProcessTransformToForegroundApplication);

	MainFrame *mFrame = new MainFrame(wxT("Active Graph Interface 3D"));
	//mFrame->initload();
	mFrame->Show(true);
	return true;
}
