#pragma once
#include <queue>
#include "CNTKLibrary.h"
#include <string>
#include <array>
#include <stack>
#include <set>
using namespace std;
using namespace CNTK;

struct nodeType {
	string s;
	int rythmIndex;
	int syllableNumber;
	int accentPos;
	int endIndex;
	int index;
	nodeType();
	nodeType(string s);
	nodeType(int p);
	string end();
	static array<array<unsigned char, 256>,2> lettersToReplace;
	static array<unsigned char, 256> vowelsToReplace;
	static std::unordered_map<string, size_t> endStringToCode;
	static vector<string> endCodeToString;
	static array<array<int, 10>, 10> rythmToIndex; // accentPos, syllNumber to index. Should be created in poem constructor
	static unordered_map<string, int> accentizer;
};

struct wordsType {
	int number;
	static std::vector<float> slotSeq;
	friend bool operator<(const wordsType& a, const wordsType& b) {
		return slotSeq[a.number] < slotSeq[b.number];
	}
};

struct linePatternType {
	vector<int> sylls;
	vector<int> weight;
	friend bool operator<(const linePatternType& a, const linePatternType& b) {
		for (int i = 0; i < a.weight.size(); i++) {
			if (a.weight[i] > b.weight[i]) {
				return true;
			}
			else if (a.weight[i] < b.weight[i])
				return false;
		}
		return false;
	}
};

struct patternType {
	vector<int> weight;
	priority_queue<linePatternType> nextRythmedPatternsPQ[4]; // separate queue for every line
	stack<linePatternType> prevRythmedPatterns[4];
	array<set<int>,2> ends;
	void nextPattern(int lineIndex);
	void prevPattern(int lineIndex);
	friend bool operator<(const patternType& a, const patternType& b) {
		for (int i = 0; i < a.weight.size(); i++) {
			if (a.weight[i] > b.weight[i]) {
				return true;
			}
			else if (a.weight[i] < b.weight[i])
				return false;
		}
		return false;
	}
};

class poem
{
public:
	poem();
	string excludeWord(int sentence, int word);
	vector<string> possibleWords(int slot, int line);
	string nextFour(string s, int rhythmLenth, int rhythmRegularity);
	string prevFour();
	static array<bool,256> vowelsRus;
	struct syllableType {
		int total;
		int regularity;
		int firstAccent;
	};
	static syllableType syllableInfo;
	static vector<array<int, 2>> rythms; // firstAccent, syllableNumber. All are unique. Should be created in poem constructor
	static std::unordered_map<string, size_t> vocabToIndex;
	static int indexOf$;
	set<string> possibleEnds(int index);
	string updateEnd(int endIndex, string end);
	string nextPattern(int lineIndex);
	string prevPattern(int lineIndex);
	array<string, 2> currentEnds();
private:
	string text;
	bool updateText(string s);

	string nextWord();
	std::priority_queue< wordsType> pqWords[2]; // queue of all words

	static FunctionPtr nextParModelFunc;
	static size_t vocabSize;
	static Variable nextParInputVarForward;
	static Variable nextParInputVarBackward;
	static int indexOfDot;
	static int indexOfSpace;
	vector<nodeType> nodes;

	size_t outputSampleSize;

	vector < map<int, map<int, set<int>>>> syllAndEndToWords; // for every slot - mapping of syll/end to words

	vector<int> rythmPattern;
	array<int,2> ends;
	int slotOfNewWord;
	bool rythmed(int rythmedBefore, int ap, int sn);
	int addedWord;
	int oppositeSlot;
	int rythmIndexOfAddedWord;
	vector<int> addedWords;
	vector<vector<int>> resultVectors;
	void print2();

	vector < map<int, map<int, vector<int>>>> syllEndToWordsOrdered; // for every slot - mapping of syll/end to words
	vector < map<int, vector<int>>> syllToWordsOrdered; // for every slot - mapping of syll/end to words
	vector<int> orderedWords;
	priority_queue<patternType> nextRythmedPatternsPQ; 
	stack<array<vector<int>,4>> prevRythmedPatterns;

	string result();
	string bestSentence();

	unordered_map<int, int> wordToOrderedIndex;

	unordered_map<int, map<int, set<int>>> endToSlotToSyll;

	map<int, map<int, map<int, map<int, set<int>>>>> pairs; // slot0 slot1 syll0 syll1 end

	map<int, map <int, map<int, map<int, array<int, 3>>>>> fours[4]; // lineIndex - slot0 - syll0 - slot1 - syll1 - end scenario (original/changed) - offsets permitted (links to starting slots for three dot scenarios)
	map<int, vector<array<map<int,int>,2>>> accumulatedLengths; // slot0-slot1-alti-consistent/non-consistent to all next words. Alti=1 starts from original or zero, alti=0 - from all possible.
	bool updateFourForNewEndRecoursively(int lineIndex, bool prevLineScenario);
	void updateAccumLengthsForNewEnd(vector<array<map<int,int>, 2>> *al, int startSlotOfALVector, int endSlot);
	void updateFourForNewSyllRecoursively(int lineIndex, int slot, bool dotInTheMiddle);
	vector <int> slots;
	vector <int> sylls;
	void findRythmedPatternsBasedOnFour();
	void updateAccumLengthsForNewSyllRecoursively(vector<array<map<int,int>, 2>> *al, int slot, int locSlot, const map<int, int>*  prevSet, int alti);

	vector<array<int,4>> scenarios;
	array<int, 4> currentScenarios;

	vector<int>nextDotPos;
	vector<int>prevDotPos;

	int endIndexOfAddedWord;
	void fillNextAccLen(map<int, int> *prev, int syll, map<int, int> *res, int limitedConsistency);
	bool accLenIsPermitted(vector<array<map<int, int>, 2>> *al, int syll, int prevAccLen, int locSlot, int accLenTableIndex, int limitedConsistency);
	map<int, int> startSet;// = { 0 };
	vector<bool> zeroPrevPrev;
	vector<bool> zeroCurr;
	vector<bool> zeroPrev;

	map<int, map<int, set<int>>> zeroMap;// = { 0,{} };
	int accLenTableIndex(int lineIndex);

	bool alternativeEndsPermitted[2];
	vector<bool> endSlotPermitted;

	vector<bool> spacePermitted;
	vector<int>nextDotPosPossible;
	//vector<int>prevDotPosPossible;

	stack<patternType> prevRythmedFours;
	string bestSentenceOfBestFour();
	void updatePermittedEnds(int slot0, int slot1, int syll1, map<int, set<int>> &permittedEnds, int slot, map<int, array<int, 3>> *scenTables);
	int wordCounter;

	string bestSentenseForPattern();
	patternType currentPattern;

	float PQWeightRatio;
	void updateRythm(int rhythmLenth, int rhythmRegularity);
};
