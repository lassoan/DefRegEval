#include "vtkKWApplication.h"
#include "vtkKWMenu.h"
#include "vtkKWWindowBase.h"
#include "vtkKWProstateErrorMapRenderingWidget.h"

//This is the entry point for the program.

//This function creates the window for the simulator and gives
//control to its event handler. The program terminates after
//the user closes the window.

int main(int argc, char *argv[])
{

	// Initialize Tcl
	Tcl_Interp *interp = vtkKWApplication::InitializeTcl(argc, argv, &cerr);
	if (!interp)
	{
		cerr << "Error: InitializeTcl failed" << endl ;
		return 1;
	}

	//create directory to save settings file
	CreateDirectory("C:\\ProstateErrorMapRendering", NULL);

	//Create the application
	vtkKWApplication *app = vtkKWApplication::New();
	app->SetName("ProstateErrorMapRendering");
	app->SetRegistryLevel(0);
    app->PromptBeforeExitOff();
	app->RestoreApplicationSettingsFromRegistry();

	//Create the window
	vtkKWWindowBase *win = vtkKWWindowBase::New();
	app->AddWindow(win);
	win->Create();
	win->SetPosition(50, 50);
	win->SetSize(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT);
	win->GetMenu()->DeleteAllItems();
	win->SetResizable(false, false);

	//add the widget
	vtkKWProstateErrorMapRenderingWidget* widget = vtkKWProstateErrorMapRenderingWidget::New();

	widget->SetParent(win->GetViewFrame());
	widget->SetParentWindow(win);
	widget->Create();

	app->Script("pack %s -expand y -fill both -anchor c -expand y", 
              widget->GetWidgetName());

	//Start the application
	int ret = 0;
	win->Display();

	app->Start(argc, argv);
	ret = app->GetExitStatus();

	win->Close();

	//Deallocate and exit
	widget->Delete();
	win->Delete();
	app->Delete();

	return ret;
}
