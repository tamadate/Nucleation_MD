# pragma once
#include "../../functions.hpp"


class StickPosition : public Functions{
public:
	void postLoop(void);
    void initial(Variables *VARS, Physical *PP, MDcondition *CON, Observer *OBS);
    std::vector<int> stickPositionList;

    char filename[100];

	StickPosition(std::vector<string> readings){
        type=0;
        step=stoi(readings[3]);
        char vaporStickFile[100];
        strcpy(vaporStickFile,readings[2].c_str());
        ifstream stream(vaporStickFile);
        string str;
        while(getline(stream,str)) {
            if(str.length()==0) continue;
            stickPositionList.push_back(stoi(str)-1);
        }
        cout<<"Vapor stick position file -->\t\t"<<readings[2]<<endl;
    };
	~StickPosition(void){};

private:

};
