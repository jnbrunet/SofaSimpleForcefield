// Here are just several convenient functions to help user to know what contains the plugin

extern "C" {
void        initExternalModule();
const char* getModuleName();
const char* getModuleVersion();
const char* getModuleLicense();
const char* getModuleDescription();
const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return "SofaSimpleForcefield";
}

const char* getModuleVersion()
{
    return "1.0";
}

const char* getModuleLicense()
{
    return "LGPL 2.1";
}

const char* getModuleDescription()
{
    return "Simple implementation of a Saint-Venant-Kirchhoff force field";
}

const char* getModuleComponentList()
{
    return "SVKElasticForcefield";
}