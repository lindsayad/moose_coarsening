#include "coarseningApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<coarseningApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

coarseningApp::coarseningApp(InputParameters parameters) : MooseApp(parameters)
{
  coarseningApp::registerAll(_factory, _action_factory, _syntax);
}

coarseningApp::~coarseningApp() {}

void
coarseningApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"coarseningApp"});
  Registry::registerActionsTo(af, {"coarseningApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
coarseningApp::registerApps()
{
  registerApp(coarseningApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
coarseningApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  coarseningApp::registerAll(f, af, s);
}
extern "C" void
coarseningApp__registerApps()
{
  coarseningApp::registerApps();
}
