#include "wx/wx.h"
#include "wx/splitter.h"
#include "wx/panel.h"
#include <wx/menu.h>
#include "DrawPanel.h"
#include "SubPanel.h"
#include "AppearanceWindow.h"

// MainFrame extends wxFrame
class MainFrame : public wxFrame{
public:
	MainFrame(const wxString& title);

	wxSplitterWindow *base;
	DrawPanel *left;
	SubPanel *right;

	wxMenuBar *menubar;
	wxMenu *file;
	wxMenu *demoMenu;
	wxMenu *demoMenu2;
	wxMenu *edit;
	wxMenu *layoutMenu;
	wxMenu *rotationMenu;
	wxMenuItem *quit;

	AppearanceWindow * appw;

	void initload();
	//void GetSubPanel(SubPanel *);
	//void GetDrawPanel(DrawPanel *);
	SubPanel * GetSubPanel();
	DrawPanel * GetDrawPanel();
	void ResetMenuParams();
	void OnQuit(wxCommandEvent& event);
	void CaptureImage(wxCommandEvent& event);
	void Reset(wxCommandEvent& event);
	void ChangeLayoutModeTo2D(wxCommandEvent& event);
	void ChangeLayoutModeTo3D(wxCommandEvent& event);
	void SetAutoXRotation(wxCommandEvent& event);
	void SetAutoYRotation(wxCommandEvent& event);
	void StopAutoRotation(wxCommandEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OpenAppearanceWindow(wxCommandEvent& event);

	void LoadDemoData1(wxCommandEvent& event);
	void LoadDemoData2(wxCommandEvent& event);
	void LoadDemoData3(wxCommandEvent& event);
	void LoadDemoData4(wxCommandEvent& event);
	void LoadDemoData5(wxCommandEvent& event);

	void LoadDemoData6(wxCommandEvent& event);
	void LoadDemoData7(wxCommandEvent& event);
	void LoadDemoData8(wxCommandEvent& event);
	void LoadDemoData9(wxCommandEvent& event);
	void LoadDemoData10(wxCommandEvent& event);
};