#include <vector>
#include <iostream>
#include <fstream>
#include <optional>
#include <string>
#include <sstream>
#include <stdexcept>
#include "../filters/include/filter.hpp"

std::vector<std::string> readFileToVector(const std::string& filename)
{
    std::vector<std::string> lines;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error cant open file " << filename << "'" << std::endl;
        return lines;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    
    file.close();
    return lines;
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::vector<std::string> lines = readFileToVector(argv[1]);
    if (lines.empty()) {
        std::cerr << "No data read from input file!" << std::endl;
        return 1;
    }

    filters::filterEMF<double> EMF(filters::EMFenvironment::PHYSICALS, 0.2, 0.0, 0.0, 0.0);
    filters::filterMedian<double> Med(16);
    std::vector<double> aXV, aYV, aZV, wXV, wYV, wZV;
    
    // ИСПРАВЛЕНО: правильное открытие файла для записи
    std::ofstream outputFile("output.txt", std::ios::out);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Cannot open output.txt for writing!" << std::endl;
        return 1;
    }

    for (const auto& line : lines)
    {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string tS, wX, wY, wZ, aX, aY, aZ;
        
        std::getline(ss, tS, ',');
        std::getline(ss, wX, ',');
        std::getline(ss, wY, ',');
        std::getline(ss, wZ, ',');
        std::getline(ss, aX, ',');
        std::getline(ss, aY, ',');
        std::getline(ss, aZ, ',');
        
        double timeStamp, aXValue, aYValue, aZValue, wXValue, wYValue, wZValue;
        try
        {
            timeStamp = std::stod(tS);
            wXValue = std::stod(wX);
            wYValue = std::stod(wY);
            wZValue = std::stod(wZ);
            aXValue = std::stod(aX);
            aYValue = std::stod(aY);
            aZValue = std::stod(aZ);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "Invalid data in line: " << line << " - " << e.what() << std::endl;
            continue;
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "Out of range error in line: " << line << " - " << e.what() << std::endl;
            continue;
        }

        aXV.push_back(aXValue);
        aYV.push_back(aYValue);
        aZV.push_back(aZValue);
        wXV.push_back(wXValue);
        wYV.push_back(wYValue);
        wZV.push_back(wZValue);

        if (aXV.size() > 128)
        {
            aXV.erase(aXV.begin());
            aYV.erase(aYV.begin());
            aZV.erase(aZV.begin());
            wXV.erase(wXV.begin());
            wYV.erase(wYV.begin());
            wZV.erase(wZV.begin());
        }

        outputFile << "$GYRACC";

        Med.setSignal(wXV);
        Med.applyFilter();
        std::vector<double> buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        Med.setSignal(wYV);
        Med.applyFilter();
        buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        Med.setSignal(wZV);
        Med.applyFilter();
        buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        Med.setSignal(aXV);
        Med.applyFilter();
        buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        Med.setSignal(aYV);
        Med.applyFilter();
        buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        Med.setSignal(aZV);
        Med.applyFilter();
        buffer = Med.getFilteredSignal();
        // EMF.setSignal(buffer);
        // EMF.applyFilter();
        // buffer = EMF.getFilteredSignal();
        if (!buffer.empty()) {
            outputFile << "," << buffer.back();
        } else {
            outputFile << ",0";
        }

        outputFile << "," << timeStamp / 1000.0 << std::endl;
    }

    outputFile.close();
    std::cout << "Data processing completed. Output written to output.txt" << std::endl;
    
    // Проверяем, создался ли файл
    std::ifstream checkFile("output.txt");
    if (checkFile.is_open()) {
        std::cout << "File output.txt successfully created!" << std::endl;
        checkFile.close();
    } else {
        std::cerr << "Warning: output.txt was not created!" << std::endl;
    }
    
    return 0;
}