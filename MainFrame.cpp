#include "MainFrame.h"
#include <fstream>
#include <iostream>
using namespace std;

void loadMatrixData_t(const char *);
void loadLayoutData_t(const char *);
void loadMatrixData_b(const char *);
void loadLayoutData_b(const char *);

extern string graphName;
static bool x_rotation;
static bool y_rotation;

MainFrame::MainFrame(const wxString& title)
: wxFrame(NULL, wxID_ANY,title,wxDefaultPosition,wxSize(1280,870)){

	base = new wxSplitterWindow(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_LIVE_UPDATE);
	base->SetSashGravity(0.80);
	base->SetMinimumPaneSize(100);

	//Panels
	int args[] = {WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16};
	left = new DrawPanel(base, args);
	right = new SubPanel(base);

	base->SplitVertically(left,right);

	//MenuBar
	menubar = new wxMenuBar();
	file = new wxMenu();
	file->Append(wxID_OPEN, wxT("Open"));
	file->AppendSeparator();

	demoMenu = new wxMenu();
	demoMenu->Append(20, wxT("lesmis"));
	demoMenu->Append(21, wxT("busho2"));
	demoMenu->Append(22, wxT("polblogs"));
	demoMenu->Append(23, wxT("math"));
	demoMenu->Append(24, wxT("ca-Astro"));
	file->AppendSubMenu(demoMenu, wxT("Demo"));

	demoMenu2 = new wxMenu();
	demoMenu2->Append(30, wxT("wiss1"));
	demoMenu2->Append(31, wxT("wiss2"));
	demoMenu2->Append(32, wxT("wiss3"));
	demoMenu2->Append(33, wxT("niconicogakkai"));
	demoMenu2->Append(34, wxT("4daigaku"));
	file->AppendSubMenu(demoMenu2, wxT("WISS Demo"));

	edit = new wxMenu();

	edit->Append(10, wxT("Reset"));
	edit->AppendSeparator();

	edit->Append(16, wxT("Appearance"));
	edit->AppendSeparator();

	layoutMenu = new wxMenu();
	layoutMenu->Append(11, wxT(" 3D  ✓"));
	layoutMenu->Append(12, wxT(" 2D"));

	edit->AppendSubMenu(layoutMenu, wxT("Layout Mode"));
	edit->AppendSeparator();

	rotationMenu = new wxMenu();
	rotationMenu->Append(13, wxT(" X Rotation  "));
	rotationMenu->Append(14, wxT(" Y Rotation  "));
	rotationMenu->Append(15, wxT(" Stop "));
	x_rotation = false;
	y_rotation = false;

	edit->AppendSubMenu(rotationMenu, wxT("Auto Rotation"));
	edit->AppendSeparator();

	edit->Append(17, wxT("Capture Image"));

	quit = new wxMenuItem(file, wxID_EXIT, wxT("&Quit\tCtrl+W"));
	file->Append(quit);

	menubar->Append(file, wxT("&File"));
	menubar->Append(edit, wxT("&Edit"));
	SetMenuBar(menubar);

	appw = new AppearanceWindow(right ,wxT("Appearance"));

	//Connect Event
	Connect(wxID_EXIT, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::OnQuit));
	Connect(10, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::Reset));
	Connect(11, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::ChangeLayoutModeTo3D));
	Connect(12, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::ChangeLayoutModeTo2D));
	Connect(13, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::SetAutoXRotation));
	Connect(14, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::SetAutoYRotation));
	Connect(15, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::StopAutoRotation));
	Connect(16, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::OpenAppearanceWindow));
	Connect(17, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::CaptureImage));

	Connect(20, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData1));
	Connect(21, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData2));
	Connect(22, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData3));
	Connect(23, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData4));
	Connect(24, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData5));

	Connect(30, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData6));
	Connect(31, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData7));
	Connect(32, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData8));
	Connect(33, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData9));
	Connect(34, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::LoadDemoData10));

	Connect(wxID_OPEN, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(MainFrame::OnOpen));

	this->Centre();
}

SubPanel * MainFrame::GetSubPanel(){
	return right;
}

DrawPanel * MainFrame::GetDrawPanel(){
	return left;
}

void MainFrame::OnQuit(wxCommandEvent& WXUNUSED(event)){
	Close(true);
}

void MainFrame::initload(){
	wxFileDialog * openFileDialog = new wxFileDialog(this);
	if (openFileDialog->ShowModal() == wxID_OK){
		wxString filePath = openFileDialog->GetPath();
		string fname = string(filePath.mb_str());
		int path_i = fname.find_last_of("/");
		int ext_i = fname.find_last_of(".");
		string f_ext = fname.substr(ext_i + 1);
		string filename = fname.substr(path_i + 1);

		if(f_ext == "txt"){
			if(filename.find("DisMat") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 11);
				loadMatrixData_t(filePath.mb_str());
			}
			else if(filename.find("Data.") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 9);
				loadLayoutData_t(filePath.mb_str());
			}
			else if(filename.find("DataAll") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 12);
				loadLayoutData_t(filePath.mb_str());
			}
			else{
				cerr << "This file is not available" << endl;
				exit(1);
			}
		}
		else if(f_ext == "bin"){
			cout << "bin" << endl;
			if(filename.find("DisMat") != string::npos){
				cout << "DisMat" << endl;
				graphName = fname.substr(path_i+1, fname.size() - path_i - 11);
				loadMatrixData_b(filePath.mb_str());
			}
			else if(filename.find("Data.") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 9);
				loadLayoutData_b(filePath.mb_str());
			}
			else if(filename.find("DataAll") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 12);
				loadLayoutData_b(filePath.mb_str());
			}
			else{
				cerr << "This file is not available" << endl;
				exit(1);
			}
		}
	    cout << "graphName:" << graphName << endl;
		left->SetupPanel();
		right->Init();
		appw->Init();
	}
}

void MainFrame::OnOpen(wxCommandEvent& event){
	wxFileDialog * openFileDialog = new wxFileDialog(this);
	if (openFileDialog->ShowModal() == wxID_OK){
		wxString filePath = openFileDialog->GetPath();
		string fname = string(filePath.mb_str());
		int path_i = fname.find_last_of("/");
		int ext_i = fname.find_last_of(".");
		string f_ext = fname.substr(ext_i + 1);
		string filename = fname.substr(path_i + 1);

		if(f_ext == "txt"){
			if(filename.find("DisMat") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 11);
				loadMatrixData_t(filePath.mb_str());
			}
			else if(filename.find("Data.") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 9);
				loadLayoutData_t(filePath.mb_str());
			}
			else if(filename.find("DataAll") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 12);
				loadLayoutData_t(filePath.mb_str());
			}
			else{
				cerr << "This file is not available" << endl;
				exit(1);
			}
		}
		else if(f_ext == "bin"){
			if(filename.find("DisMat") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 11);
				loadMatrixData_b(filePath.mb_str());
			}
			else if(filename.find("Data.") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 15);
				loadLayoutData_b(filePath.mb_str());
			}
			else if(filename.find("DataAll") != string::npos){
				graphName = fname.substr(path_i+1, fname.size() - path_i - 18);
				loadLayoutData_b(filePath.mb_str());
			}
			else{
				cerr << "This file is not available" << endl;
				exit(1);
			}
		}
    	cout << "graphName:" << graphName << endl;
		left->SetupPanel();
		right->Init();
		appw->Init();
	}
}

void MainFrame::ResetMenuParams(){
	x_rotation = false;
	left->SetXRotation(false);
	rotationMenu->SetLabel(13, wxT(" X Rotation"));
	y_rotation = false;
	left->SetYRotation(false);
	rotationMenu->SetLabel(14, wxT(" Y Rotation"));
}

void MainFrame::Reset(wxCommandEvent& event){
	left->ResetLayout();
	right->Init();
	appw->Init();
	ResetMenuParams();
}

void MainFrame::ChangeLayoutModeTo2D(wxCommandEvent& event){
	left->ChangeLayoutMode(2);
	right->Init();
	appw->Init();
	layoutMenu->SetLabel(11, wxT(" 3D"));
	layoutMenu->SetLabel(12, wxT(" 2D  ✓"));
	ResetMenuParams();
}

void MainFrame::ChangeLayoutModeTo3D(wxCommandEvent& event){
	left->ChangeLayoutMode(3);
	right->Init();
	appw->Init();
	layoutMenu->SetLabel(11, wxT(" 3D  ✓"));
	layoutMenu->SetLabel(12, wxT(" 2D"));
	ResetMenuParams();
}

void MainFrame::SetAutoXRotation(wxCommandEvent& event){
	if(x_rotation){
		left->SetXRotation(false);
		x_rotation = false;
		rotationMenu->SetLabel(13, wxT(" X Rotation  "));
	}
	else{
		left->SetXRotation(true);
		x_rotation = true;
		rotationMenu->SetLabel(13, wxT(" X Rotation  ✓"));
	}
}

void MainFrame::SetAutoYRotation(wxCommandEvent& event){
	if(y_rotation){
		left->SetYRotation(false);
		y_rotation = false;
		rotationMenu->SetLabel(14, wxT(" Y Rotation  "));
	}
	else{
		left->SetYRotation(true);
		y_rotation = true;
		rotationMenu->SetLabel(14, wxT(" Y Rotation  ✓"));
	}
}

void MainFrame::StopAutoRotation(wxCommandEvent& event){
	ResetMenuParams();
}

void MainFrame::CaptureImage(wxCommandEvent& event){
	left->SavePixelData();
}

void MainFrame::OpenAppearanceWindow(wxCommandEvent& event){
	appw->Show(true);
	appw->Raise();
}

void MainFrame::LoadDemoData1(wxCommandEvent& event){
	graphName = "lesmis";
	string filePath = "../data/lesmis/lesmisDisMat.txt";
	loadMatrixData_t(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData2(wxCommandEvent& event){
	graphName = "busho2";
	string filePath = "../data/busho2/busho2DisMat.txt";
	loadMatrixData_t(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData3(wxCommandEvent& event){
	graphName = "polblogs";
	string filePath = "../data/polblogs/polblogsLayoutDataAll.bin";
	loadLayoutData_b(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData4(wxCommandEvent& event){
	graphName = "math";
	string filePath = "../data/math/mathLayoutDataAll.bin";
	loadLayoutData_b(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData5(wxCommandEvent& event){
	graphName = "ca-Astro";
	string filePath = "../data/ca-Astro/ca-AstroLayoutDataAll.bin";
	loadLayoutData_b(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData6(wxCommandEvent& event){
	graphName = "wiss";
	string filePath = "../data/wiss/wissDisMat.txt";
	loadMatrixData_t(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData7(wxCommandEvent& event){
	graphName = "wiss2";
	string filePath = "../data/wiss2/wiss2DisMat.txt";
	loadMatrixData_t(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData8(wxCommandEvent& event){
	graphName = "wiss3";
	string filePath = "../data/wiss3/wiss3DisMat.txt";
	loadMatrixData_t(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData9(wxCommandEvent& event){
	graphName = "niconicogakkai";
	string filePath = "../data/niconicogakkai/niconicogakkaiLayoutDataAll.bin";
	loadLayoutData_b(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}

void MainFrame::LoadDemoData10(wxCommandEvent& event){
	graphName = "4daigaku";
	string filePath = "../data/4daigaku/4daigakuLayoutDataAll.bin";
	loadLayoutData_b(filePath.c_str());
	left->SetupPanel();
	right->Init();
	appw->Init();
}