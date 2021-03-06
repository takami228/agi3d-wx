#include "SubPanel.h"
#include "MainFrame.h"
#include <iostream>
#include <string>
#include <vector>
#include <cml/cml.h>
#include <map>

typedef cml::vector4f vector4;
using namespace std;

static DrawPanel* dp;
static float fps = 0.0f;
static float tmp = 0.0f;
static wxStaticText * m_FPS;
static wxStaticText * fileNameLabel;
static wxStaticText * nodeSizeLabel;
static wxStaticText * edgeSizeLabel;

static int nodeThresholdAttrID = 0;
static int edgeThresholdAttrID = 0;

static wxListBox * target;
static wxListBox * listbox;
static int pickid = -1;

static float _nodethreshold_t = 10000000;
static float _nodethreshold_b = -1;
static map<int,int> labelMap;
static int node_slider_b_pos = 0;
static int node_slider_t_pos = 100;
static int edge_slider_b_pos = 0;
static int edge_slider_t_pos = 100;

extern vector<int> * neighbor;
extern vector<string> labels;
extern float * nodevalues;
extern float nodevalue_max;
extern float nodevalue_min;
extern string filename;
extern int N;
extern int M;
extern string graphName;

static int _id = -1;

void loadNodeAttrData(int);
void loadEdgeAttrData(int);

BEGIN_EVENT_TABLE(SubPanel, wxPanel)
EVT_LISTBOX(555, SubPanel::handleListEvent)
END_EVENT_TABLE()

SubPanel::SubPanel(wxWindow* parent)
: wxPanel(parent, wxID_ANY, wxDefaultPosition, wxSize(280,870)){

	wxWindow * p =  this->GetParent();

	MainFrame * mf = (MainFrame *)(p->GetParent());
	dp = (DrawPanel *)mf->GetDrawPanel();

	wxPanel *myPanel = new wxPanel(this, wxID_ANY, wxDefaultPosition, wxSize(280,870));

	wxBoxSizer *mybox = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer *b1 = new wxBoxSizer(wxHORIZONTAL);
	//wxImage::AddHandler( new wxPNGHandler );

	m_FPS = new wxStaticText(this, -1, wxT("   FPS"), wxDefaultPosition, wxSize(100,28), wxALIGN_LEFT);
	fileNameLabel = new wxStaticText(this, -1, wxT("   filename"), wxDefaultPosition, wxSize(100,25), wxALIGN_LEFT);
	nodeSizeLabel = new wxStaticText(this, -1, wxT("   #Node "), wxDefaultPosition, wxSize(100,25),wxALIGN_LEFT);
	edgeSizeLabel = new wxStaticText(this, -1, wxT("   #Edge  "), wxDefaultPosition, wxSize(100,25),wxALIGN_LEFT);

	mybox->Add(m_FPS, 1, wxEXPAND, 5);
	mybox->Add(fileNameLabel, 1, wxEXPAND, 0);
	mybox->Add(nodeSizeLabel, 1, wxEXPAND, 0);
	mybox->Add(edgeSizeLabel, 1, wxEXPAND, 0);

	mybox->Add(myPanel, 1, wxEXPAND | wxALL,5);
	this->SetSizer(mybox);

	/*
	wxStaticText * cameraPositionText = new wxStaticText(myPanel, wxID_ANY, wxT("Zoom In/Out"), wxPoint(10, 0));
	cameraPositionSlider = new wxSlider(myPanel, 30, 100, 10, 200, wxPoint(10, 20), wxSize(165, -1));
	*/

	wxStaticText * DeltaText = new wxStaticText(myPanel, wxID_ANY, wxT("Projection Factor"), wxPoint(10, 0));
  	DeltaSlider = new wxSlider(myPanel, 60, 50, 0, 400, wxPoint(10, 20), wxSize(245, -1));

	wxStaticText * nodeThresholdText = new wxStaticText(myPanel, wxID_ANY, wxT("Node Filtering"), wxPoint(10,50));

	wxArrayString nodeAttrs;
	nodeAttrs.Add("No Feature");
	nodeAttrs.Add("Degree Centrality");
	nodeAttrs.Add("Closness Centrality");
	nodeAttrs.Add("Betweeness Centrality");
	nodeAttrs.Add("Clustering Coefficient");
	nodeAttrs.Add("Eigenvector Centrality");
	nodeAttrs.Add("PageRank");
	nodeAttrsChoice = new wxChoice(myPanel, 120, wxPoint(7, 70), wxSize(250,-1), nodeAttrs);

	nodeThresholdSlider_b = new wxSlider(myPanel, 70, 0, 0, 100, wxPoint(10, 95), wxSize(247, -1));
	nodeThresholdSlider_t = new wxSlider(myPanel, 71, 100, 0, 100, wxPoint(10, 112), wxSize(247, -1));
	node_slider_b_pos = 0;
	node_slider_t_pos = 100;

	wxStaticText * edgeThresholdText = new wxStaticText(myPanel, wxID_ANY, wxT("Edge Filtering"), wxPoint(10, 140));

	wxArrayString edgeAttrs;
	edgeAttrs.Add("No Feature");
	edgeAttrs.Add("Simpson");
	edgeAttrs.Add("Extended Simpson");
	edgeAttrs.Add("Betweeness");
	edgeAttrs.Add("Weight");
	edgeAttrsChoice = new wxChoice(myPanel, 121, wxPoint(7,165), wxSize(250,-1), edgeAttrs);

	edgeThresholdSlider_b = new wxSlider(myPanel, 80, 0, 0, 100, wxPoint(10, 190), wxSize(247, -1));
	edgeThresholdSlider_t = new wxSlider(myPanel, 81, 100, 0, 100, wxPoint(10, 207), wxSize(247, -1));
	edge_slider_b_pos = 0;
	edge_slider_t_pos = 100;

	target = new wxListBox(myPanel, 554, wxPoint(5, 240), wxSize(260,25));

	listbox = new wxListBox(myPanel, 555, wxPoint(5, 270), wxSize(260,460));
	labelMap.clear();

	//Connect(30, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::UpdateCameraPosition));
	Connect(40, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::NortifyUpdateNodeSize));
	Connect(50, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::NortifyUpdateEdgeThickness));
	Connect(60, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::NortifyUpdateDelta));
	Connect(61, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::NortifyUpdateScale));
	Connect(62, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::NortifyUpdateDimension));

	Connect(70, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::UpdateNodeValueThreshold_b));
	Connect(71, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::UpdateNodeValueThreshold_t));
	Connect(80, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::UpdateEdgeValueThreshold_b));
	Connect(81, wxEVT_COMMAND_SLIDER_UPDATED, wxScrollEventHandler(SubPanel::UpdateEdgeValueThreshold_t));

	Connect(120, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(SubPanel::SelectNodeAttr));
	Connect(121, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(SubPanel::SelectEdgeAttr));

	Connect(200, wxEVT_COMMAND_CHECKBOX_CLICKED,  wxCommandEventHandler(SubPanel::OnToggleEdge));
	Connect(201, wxEVT_COMMAND_CHECKBOX_CLICKED,  wxCommandEventHandler(SubPanel::OnToggleNodeSize));

	Connect(202, wxEVT_COMMAND_CHECKBOX_CLICKED,  wxCommandEventHandler(SubPanel::OnToggleXRotation));
	Connect(203, wxEVT_COMMAND_CHECKBOX_CLICKED,  wxCommandEventHandler(SubPanel::OnToggleYRotation));
}

void SubPanel::SetFPS(float _f){
	fps = _f;
	wxString str = wxString::Format(wxT("   FPS : %f"), fps);
	m_FPS->SetLabel(str);
	tmp = _f;
}

void SubPanel::Init(){
	//cameraPositionSlider->SetValue(100);
	DeltaSlider->SetValue(50);

	nodeThresholdSlider_b->SetValue(0);
	nodeThresholdSlider_t->SetValue(100);
	float _nodethreshold_b = -1;
	float _nodethreshold_t = 10000000;

	edgeThresholdSlider_b->SetValue(0);
	edgeThresholdSlider_t->SetValue(100);

	nodeAttrsChoice->SetSelection(0);
	nodeThresholdAttrID = 0;
	edgeAttrsChoice->SetSelection(0);
	edgeThresholdAttrID = 0;

	node_slider_b_pos = 0;
	node_slider_t_pos = 100;
	edge_slider_b_pos = 0;
	edge_slider_t_pos = 100;

	target->SetLabel(wxString());
	listbox->Clear();

	labelMap.clear();

	dp->UpdateSize(1.0f);
	dp->UpdateThickness(1.0f);

	dp->SetXRotation(false);
	dp->SetYRotation(false);

	wxString _graphlabel = wxString::Format(wxT("   FileName \t %s"),graphName);
	wxString _nodeSize = wxString::Format(wxT("   #Node \t %i"),N);
	wxString _edgeSize = wxString::Format(wxT("   #Edge   \t %i"),M);

	fileNameLabel->SetLabel(_graphlabel);
	nodeSizeLabel->SetLabel(_nodeSize);
	edgeSizeLabel->SetLabel(_edgeSize);
}

void SubPanel::setTarget(int id){
	target->Clear();
	listbox->Clear();
	labelMap.clear();
	if(id != -1){
		target->Append(wxString(labels[id]));
		_id = id;
		int key = 0;
		for(int i = 0; i < neighbor[id].size(); i++){
			int _nei = neighbor[id][i];
			if(nodevalues[_nei] >= _nodethreshold_b && nodevalues[_nei] <= _nodethreshold_t){
				listbox->Append(wxString(labels[_nei]));
				labelMap.insert(map<int,int>::value_type(key,i));
				key++;
			}
		}
		listbox->SetFirstItem(0);
	}
	else{
		_id = id;
	}
}

void SubPanel::handleListEvent(wxCommandEvent& event){
	int m = event.GetInt();
	cout << _id <<  " " << labelMap[m] << endl;
	dp->changeColor(labelMap[m]);
}

void SubPanel::UpdateCameraPosition(wxScrollEvent& event){
	float v = (float)(cameraPositionSlider->GetValue())/100.0f;
	dp->UpdateEye(v);
}

void SubPanel::NortifyUpdateNodeSize(wxScrollEvent& event){
	float rate = (float)(event.GetInt())/100.0f;
	dp->UpdateSize(rate);
}

void SubPanel::NortifyUpdateEdgeThickness(wxScrollEvent& event){
	float rate = (float)(event.GetInt())/50.0f;
	dp->UpdateThickness(rate);
}

void SubPanel::NortifyUpdateDelta(wxScrollEvent& event){
	//float rate = (float)(event.GetInt())/100.0f;
	float rate = (float)(DeltaSlider->GetValue())/100.0f;
	dp->ModifyDelta(rate);
}

void SubPanel::NortifyUpdateScale(wxScrollEvent& event){
	float rate = (float)(event.GetInt())/20.0f;
	dp->ScaleLayout(rate);
}

void SubPanel::NortifyUpdateDimension(wxScrollEvent& event){
	float rate = (float)(event.GetInt())/1000.0f;
	dp->ChangeDimension(rate);
}

void SubPanel::UpdateNodeValueThreshold_b(wxScrollEvent& event){
	node_slider_b_pos = nodeThresholdSlider_b->GetValue();
	if(node_slider_b_pos <= node_slider_t_pos){
		float z = (float)(node_slider_b_pos)/100.0f;
		_nodethreshold_b = dp->UpdateNodeThreshold_b(z, nodeThresholdAttrID);
	}
	else{
		nodeThresholdSlider_b->SetValue(node_slider_t_pos-1);
		_nodethreshold_b = _nodethreshold_t*0.99;
	}
}

void SubPanel::UpdateNodeValueThreshold_t(wxScrollEvent& event){
	node_slider_t_pos = nodeThresholdSlider_t->GetValue();
	if(node_slider_b_pos <= node_slider_t_pos){
		float z = (float)(node_slider_t_pos)/100.0f;
		_nodethreshold_t = dp->UpdateNodeThreshold_t(z, nodeThresholdAttrID);
	}
	else{
		nodeThresholdSlider_t->SetValue(node_slider_b_pos+1);
		_nodethreshold_t = _nodethreshold_b*1.01;
	}
}


void SubPanel::UpdateEdgeValueThreshold_b(wxScrollEvent& event){
	edge_slider_b_pos = edgeThresholdSlider_b->GetValue();
	if(edge_slider_b_pos < edge_slider_t_pos){
		float z = (float)(edge_slider_b_pos)/100.0f;
		dp->UpdateEdgeThreshold_b(z, edgeThresholdAttrID);
	}
	else{
		edgeThresholdSlider_b->SetValue(edge_slider_t_pos-1);
	}
}

void SubPanel::UpdateEdgeValueThreshold_t(wxScrollEvent& event){
	edge_slider_t_pos = edgeThresholdSlider_t->GetValue();
	if(edge_slider_b_pos < edge_slider_t_pos){
		float z = (float)(edge_slider_t_pos)/100.0f;
		dp->UpdateEdgeThreshold_t(z, edgeThresholdAttrID);
	}
	else{
		edgeThresholdSlider_t->SetValue(edge_slider_b_pos+1);
	}
}

void SubPanel::OnToggleEdge(wxCommandEvent& event){
	dp->DrawEdge();
}

void SubPanel::OnToggleNodeSize(wxCommandEvent& event){
	dp->NodeModeChange();
}

void SubPanel::OnToggleXRotation(wxCommandEvent& event){
	dp->SetXRotation(true);
}

void SubPanel::OnToggleYRotation(wxCommandEvent& event){
	dp->SetYRotation(true);
}

void SubPanel::SelectNodeAttr(wxCommandEvent& event){
	wxString label = nodeAttrsChoice->GetStringSelection();
	int n = nodeAttrsChoice->GetSelection();
	nodeThresholdSlider_b->SetValue(0);
	nodeThresholdSlider_t->SetValue(100);
	node_slider_b_pos = 0;
	node_slider_t_pos = 100;
	nodeThresholdAttrID = n;
	string strs = string(label.mb_str());
	loadNodeAttrData(n);
	dp->ResetIsDrawingNodes();
}

void SubPanel::SelectEdgeAttr(wxCommandEvent& event){
	wxString label = edgeAttrsChoice->GetStringSelection();
	int n = edgeAttrsChoice->GetSelection();
	edgeThresholdSlider_b->SetValue(0);
	edgeThresholdSlider_t->SetValue(100);
	edge_slider_b_pos = 0;
	edge_slider_t_pos = 100;
	edgeThresholdAttrID = n;
	string strs = string(label.mb_str());
	cout << n << " " << strs << endl;
	loadEdgeAttrData(n);
	dp->ResetIsDrawingEdges();
}