#include <iostream>
#include <fstream>
#include <stdexcept>

#include <cxxopts.hpp>

#include <SpiceUsr.h>

using namespace std;


int main(int argc, char** argv) {
    string exename(argv[0]);
    int verbose = 0;

    string inputFile;

    cxxopts::Options options(argv[0], " - Star Database preprocessor");
    options.add_options()
        ("v,verbose", "Verbose output", cxxopts::value<int>(verbose)->default_value("0")->implicit_value("1"))
        ("h,help", "Print help")
        ("input", "Input file (Spice bdb file)", cxxopts::value<std::string>())
        ("o,output", "Output file", cxxopts::value<std::string>())
        ("l,list","List star details of star nr #", cxxopts::value<int>())
        ("m,magnitude", "Display stars with higher mag than given", cxxopts::value<double>())
        ("positional",
          "Positional arguments: input file, output file", cxxopts::value<std::vector<std::string>>())
        ;

    options.parse_positional({"input", "output", "positional"});

    options.parse(argc, argv);

    if(options.count("help"))
    {
        cout << options.help({""}) << endl;
        return 0;
    }

    if(options.count("verbose"))
        verbose = options["verbose"].as<int>();

    if(verbose>0)
        cout << "Verbose mode: " << verbose << endl;

    if (options.count("input"))
    {
        inputFile = options["input"].as<std::string>();
    } else
    {
        cout << options.help({""}) << endl;
        return -1;
    }
        

    string ofile;
    bool write_ofile = false;
    if (options.count("output"))
    {
        ofile = options["output"].as<std::string>();
        std::cout << "Output = " << options["output"].as<std::string>() << std::endl;
        write_ofile = true;
    }

    if (options.count("positional"))
    {
        std::cout << "Positional = {";
        auto& v = options["positional"].as<std::vector<std::string>>();
        for (const auto& s : v) {
             std::cout << s << ", ";
        }
        std::cout << "}" << std::endl;
    }


    // Try to read the input file
    // Needs to be a Spice type 1 master star database
    // (EK file)
    if(verbose)
        cout << "Opening " << inputFile << " for reading.." << endl;
    //furnsh_c(inputFile.c_str()); // With defaiul settings, Spice will fail here if there is somthing wrong with the input file
    
    // Lets open the file for record oriented (low level) reading:
    SpiceInt handle = 0;
    ekopr_c(inputFile.c_str(), &handle);
    SpiceInt nseg = eknseg_c(handle);
    std::vector<int> rows;
    if(verbose>0)
        cout << inputFile << " has " << nseg << " segments." << endl;
   
    int numStars = 0;

    // Parse, check (and summarize if v>0) segments
    // Gte total number of stars
    for(auto i = 0; i < nseg; ++i)
    {
        SpiceEKSegSum ss;
        ekssum_c(handle, i, &ss);

        numStars += ss.nrows;
        rows.push_back(ss.nrows);

        if(verbose>1)
        {
            cout << "==========================" << endl;
            cout << "Segment " << i+1 << endl;
            cout << "Table name: " << ss.tabnam << endl;
            cout << "Num rows:   " << ss.nrows << endl;
            cout << "Num cols:   " << ss.ncols << endl;
            cout << endl;
        }
         
        if(ss.ncols!=15)
        {
            cout << "Error; improperly formatted star database type 1. Expected 15 columns. Try with --verbose for more information." << endl;
            return -1;
        }

        for(auto j = 0; j < ss.ncols; ++j)
        {
            if(verbose>2)
            {
                cout << "    Column:     " << ss.cnames[j] << endl;
                cout << "    Data type:  " << ss.cdescrs[j].dtype << endl;
                if ( ss.cdescrs[j].size >= 0 )
                {
                    cout << "    Dimension:  " << ss.cdescrs[i].size << endl;
                } else
                {
                    cout << "    Dimension:  Variable" << endl;
                }
                                         
                if( ss.cdescrs[j].dtype == SPICE_CHR )
                {
                    if(ss.cdescrs[j].strlen>=0)
                        cout << "    String length: " << ss.cdescrs[j].strlen << endl;
                    else
                        cout << "    String lenght: Variable " << endl;
                }

                if( ss.cdescrs[j].indexd )
                    cout << "    Indexed" << endl;

                if( ss.cdescrs[j].nullok )
                    cout << "    NULLs allowed " << endl;

                cout << endl;
            }

        }

        
    }
    
    if(verbose>0)
    {
        cout << "Found " << numStars << " stars " << endl;
        cout << "Reading star data..." << endl;
    }

    // dimension vectors for the star database size:
    std::vector<long> catalog_number;
    std::vector<double> ra;
    std::vector<double> dec;
    std::vector<double> ra_epoch;
    std::vector<double> dec_epoch;
    std::vector<double> ra_sigma;
    std::vector<double> dec_sigma;
    std::vector<double> ra_pm;
    std::vector<double> dec_pm;
    std::vector<double> ra_pm_sigma;
    std::vector<double> dec_pm_sigma;
    std::vector<double> parlax;
    std::vector<double> visual_magnitude;
    std::vector<long> dm_number;
    std::vector<std::string> spectral_type;

    // read data
    for(auto segno = 0; segno < nseg; ++segno)
    {
        // Temporaries
        SpiceDouble dval;
        SpiceInt ival;
        SpiceInt nvals;
        SpiceBoolean isnull;
        SpiceChar cvals[5];

        for(auto row = 0; row < rows[segno]; ++row)
        {
            ekrcei_c(handle, segno, row, "CATALOG_NUMBER",
                &nvals, &ival, &isnull);
            catalog_number.push_back((long)ival);
        
            ekrced_c(handle, segno, row, "RA",
                &nvals, &dval, &isnull);
            ra.push_back(dval);
            ekrced_c(handle, segno, row, "DEC",
                &nvals, &dval, &isnull);
            dec.push_back(dval);
            ekrced_c(handle, segno, row, "RA_EPOCH",
                &nvals, &dval, &isnull);
            ra_epoch.push_back(dval);
            ekrced_c(handle, segno, row, "DEC_EPOCH",
                &nvals, &dval, &isnull);
            dec_epoch.push_back(dval);
            ekrced_c(handle, segno, row, "RA_SIGMA",
                &nvals, &dval, &isnull);
            ra_sigma.push_back(dval);
            ekrced_c(handle, segno, row, "DEC_SIGMA",
                &nvals, &dval, &isnull);
            dec_sigma.push_back(dval);
            ekrced_c(handle, segno, row, "RA_PM",
                &nvals, &dval, &isnull);
            ra_pm.push_back(dval);
            ekrced_c(handle, segno, row, "DEC_PM",
                &nvals, &dval, &isnull);
            dec_pm.push_back(dval);
            ekrced_c(handle, segno, row, "RA_PM_SIGMA",
                &nvals, &dval, &isnull);
            ra_pm_sigma.push_back(dval);
            ekrced_c(handle, segno, row, "DEC_PM_SIGMA",
                &nvals, &dval, &isnull);
            dec_pm_sigma.push_back(dval);
         // Parallax         double
            ekrced_c(handle, segno, row, "PARLAX",
                &nvals, &dval, &isnull);
            parlax.push_back(dval);
         // DM number        long (signed, 8 bytes)
            ekrcei_c(handle, segno, row, "DM_NUMBER",
                &nvals, &ival, &isnull);
            dm_number.push_back((long)ival);
         // Visual Magnitude double
            ekrced_c(handle, segno, row, "VISUAL_MAGNITUDE",
                &nvals, &dval, &isnull);
            visual_magnitude.push_back(dval);
         // Spectral type    char * 4
            ekrcec_c(handle, segno, row, "SPECTRAL_TYPE",
                5, &nvals, cvals, &isnull);
            spectral_type.push_back(std::string(cvals));                      
        }
    }


    if(write_ofile)
    {
        // Make sure SpiceDouble == double
        if(sizeof(SpiceDouble) != sizeof(double))
            throw runtime_error("SpiceDouble not equal to double");

        if(verbose>0)
            cout << "Writing star database to " << ofile << endl;

        ofstream fo;
        fo.open(ofile.c_str(), ios::out | ios::binary);


        // For documentation of the star data base file format,
        // see Stars.h in rtsim.

        // Write magic numbers for version 10000,10001,3
        int magicNumbers[3] = {10000, 10001, 3};
        int dummy;
        fo.write((char*)magicNumbers, 3*sizeof(int));
        
        // Star db type(1)
        dummy = 1;
        fo.write((char*)&dummy, sizeof(int));

        // Write catname lenght and catname
        // Take the cat name from the first segment
        // TODO:
        // Make guard against longer names than 127 characters
        SpiceEKSegSum ss;
        ekssum_c(handle, 0, &ss);
        string catname(ss.tabnam);
        dummy = catname.length();
        cout << "\"" << catname << "\"(" << catname.length() << ")" << endl;
        fo.write((char*)&dummy, sizeof(int));
        fo.write(catname.c_str(), dummy);

        // write magic number (10000)
        fo.write((char*)magicNumbers, sizeof(int));

        // Number of stars:
        dummy = numStars;
        fo.write((char*)&dummy, sizeof(int));

        // according to the stardb, we now write numStars number of entries
        // of each parameter, in the following order:
        // Parameter        Type
        // catalogue number long (signed, 8 bytes)
        // RA base          double
        // Dec base         double
        // RA epoch         double
        // DEC epoch        double
        // RA sigma         double
        // DEC sigma        double
        // RA pm            double
        // DEC pm           double
        // RA pm sigma      double
        // DEC pm sigma     double
        // Parallax         double
        // DM number        long (signed, 8 bytes)
        // Visual Magnitude double
        // Spectral type    char * 4

        // Write data: 
        // catalogue number long (signed, 8 bytes)
        fo.write((char*)catalog_number.data(),sizeof(long)*numStars);
        fo.write((char*)ra.data(),sizeof(double)*numStars);
        fo.write((char*)dec.data(),sizeof(double)*numStars);
        fo.write((char*)ra_epoch.data(),sizeof(double)*numStars);
        fo.write((char*)dec_epoch.data(),sizeof(double)*numStars);
        fo.write((char*)ra_sigma.data(),sizeof(double)*numStars);
        fo.write((char*)dec_sigma.data(),sizeof(double)*numStars);
        fo.write((char*)ra_pm.data(),sizeof(double)*numStars);
        fo.write((char*)dec_pm.data(),sizeof(double)*numStars);
        fo.write((char*)ra_pm_sigma.data(),sizeof(double)*numStars);
        fo.write((char*)dec_pm_sigma.data(),sizeof(double)*numStars);
        fo.write((char*)parlax.data(),sizeof(double)*numStars);
        fo.write((char*)dm_number.data(),sizeof(long)*numStars);
        fo.write((char*)visual_magnitude.data(),sizeof(double)*numStars);
        std::vector<std::string>::iterator it;
        for(it = spectral_type.begin(); it != spectral_type.end(); ++it)
        {
            fo.write((*it).c_str(),4);
        }
        fo.close();
    }
    
    if(options.count("list"))
    {
        int starno = options["list"].as<int>();
        if(starno<0 || starno>= numStars)
        {
            cout << "Error: List: Star number must be in range (0, " << numStars-1 << ")" << endl;
            return -1;
        }

        cout << "Star no\tCatnum\tRA\tDEC\tRAepoch\tDECepoch\tRAsigma\tDECsigma\tRApm\tDECpm\tRApms\tDECpms\tVisMag\tParallax\tDM_num\tSpectral" << endl;
        cout << starno << "\t";
        cout << catalog_number[starno] << "\t";
        cout << ra[starno] << "\t";
        cout << dec[starno] << "\t";
        cout << ra_epoch[starno] << "\t";
        cout << dec_epoch[starno] << "\t";
        cout << ra_sigma[starno] << "\t";
        cout << dec_sigma[starno] << "\t";
        cout << ra_pm[starno] << "\t";
        cout << dec_pm[starno] << "\t";
        cout << ra_pm_sigma[starno] << "\t";
        cout << dec_pm_sigma[starno] << "\t";
        cout << visual_magnitude[starno] << "\t";
        cout << parlax[starno] << "\t";
        cout << dm_number[starno] << "\t\"";
        cout << spectral_type[starno] << "\"" << endl; 
    }
    
    if(options.count("magnitude"))
    {
        double maglimit = options["magnitude"].as<double>();
        std::vector<double>::iterator it;
        std::vector<int> who;
        int i = 0;
        for(it = visual_magnitude.begin(); it!=visual_magnitude.end(); ++it, ++i)
            if((*it)<maglimit)
                who.push_back(i);


        cout << "Maglimit: " << maglimit << ", number of stars: " << who.size() << endl;
        if(who.size() <= 10)
        {
            cout << "Star no\tcat num\tapp mag" << endl;
            std::vector<int>::iterator jt;
            for(jt = who.begin(); jt != who.end(); ++jt)
                cout << (*jt) << "\t" << catalog_number[*jt] << "\t" << visual_magnitude[*jt] << endl;
        }
    }
                             
    return 0;
}



