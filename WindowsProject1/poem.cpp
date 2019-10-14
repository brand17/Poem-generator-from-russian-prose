#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS // "secure" CRT not available on all platforms.
#endif

#include "poem.h"
#include <sstream>
#include <fstream>
#include <algorithm>
std::vector<nodeType> indexToSlots;
const array<int, 4> slot1s = { 0,1,0,1 };
const array<int, 4> slot2s = { 2,3,2,3 };
const array<int, 4> oppositeSlot1s = { 1,0,1,0 };
const array<int, 4> oppositeSlot2s = { 3,2,3,2 };
const array<int, 4> oppositeIndices = { 2,3,0,1 };

array<bool, 256> poem::vowelsRus;
float* network1Output;
struct wordWeightType {
	int word;
	float weight;
};
vector<vector<wordWeightType>> possibleWords;
array<array<unsigned char, 256>,2> nodeType::lettersToReplace;
array<unsigned char, 256> nodeType::vowelsToReplace;
std::unordered_map<string, size_t> nodeType::endStringToCode;
vector<string> nodeType::endCodeToString;
poem::syllableType poem::syllableInfo;
FunctionPtr poem::nextParModelFunc;
//FunctionPtr poem::insertSlotModelFunc;
//FunctionPtr poem::replaceSlotModelFunc;
//FunctionPtr poem::deleteSlotModelFunc;
std::unordered_map<string, size_t> poem::vocabToIndex;
size_t poem::vocabSize;
Variable poem::nextParInputVarForward;
Variable poem::nextParInputVarBackward;
//Variable poem::insertSlotInputVarForward;
//Variable poem::insertSlotInputVarBackward;
//Variable poem::replaceSlotInputVarForward;
//Variable poem::deleteSlotInputVarForward;
int poem::indexOf$;
int poem::indexOfDot;
int poem::indexOfSpace;
array<array<int, 10>, 10> nodeType::rythmToIndex;
vector<array<int, 2>> poem::rythms;
std::vector<float> wordsType::slotSeq;
unordered_map<string, int> nodeType::accentizer;
int counter = 0;

std::shared_ptr<std::fstream> GetIfstream(const wchar_t *filePath)
{
	const size_t pathBufferLen = 1024;
	char pathBuffer[pathBufferLen];
	size_t writtenBytes = ::wcstombs(pathBuffer, filePath, pathBufferLen);
	if (writtenBytes == (size_t)-1)
		throw ("Unknown characters in the file path.");
	else if (writtenBytes == pathBufferLen)
		throw("The file path is too long");
	return std::make_shared<std::fstream>(pathBuffer);
}

nodeType::nodeType(string sss) { // to fill dictionary
	s = sss;
	syllableNumber = 0;
	int accentPosFull = 0;
	for (unsigned char i : s) {
		if (i == '`')
			accentPosFull = syllableNumber;
		else
			syllableNumber += poem::vowelsRus[i];
	}
	accentPos = accentPosFull % poem::syllableInfo.regularity;
	// find rythm index
	if (rythmToIndex[accentPos][syllableNumber] != -1) {
		rythmIndex = rythmToIndex[accentPos][syllableNumber];
	}
	else {
		rythmIndex = poem::rythms.size();
		poem::rythms.push_back({ accentPos,syllableNumber });
		rythmToIndex[accentPos][syllableNumber] = rythmIndex;
	}
	// fill end code using dictionary of ends
	string ss;// = end();
	if (syllableNumber == 0 || syllableNumber != accentPosFull+1//(syllableNumber > 3 && accentPosFull<syllableNumber-2) ||
		//0 != (syllableNumber-1-accentPos) % poem::syllableInfo.regularity
		) // check that the last vowel can be accentized
		ss = "";
	else 
		ss = end();
	auto it = endStringToCode.find(ss);
	if (it == endStringToCode.end()) {
		endStringToCode[ss] = endStringToCode.size();
		endCodeToString.push_back(ss);
	}
	endIndex = endStringToCode[ss];
	auto it2 = poem::vocabToIndex.find(sss);
	if (it2 == poem::vocabToIndex.end()) {
		index = poem::indexOf$;
	}
	else {
		index = it2->second;
	}
	// add to accent dictionary
	if (syllableNumber > 1) {
		string notAccented;
		copy_if(sss.begin(), sss.end(), back_inserter(notAccented), [](char c) {return c != '`'; });
		if(accentizer.find(notAccented)==accentizer.end()) // to ensure the best option is choosen for accentizer
			accentizer[notAccented] = index;
	}
}

nodeType::nodeType(int p) { // to initialize based on dictionary
	syllableNumber = indexToSlots[p].syllableNumber;
	accentPos = indexToSlots[p].accentPos;
}

nodeType::nodeType() {
}

string nodeType::end()
{
	std::string output;
	output.reserve(s.size()); // optional, avoids buffer reallocations in the loop
	for (size_t j = 0; j < s.size(); ++j)
		if (s.at(j) != '`')
			output += s.at(j);
	int cycles = 1;
	int last = output.length() - 1;
	unsigned char c = output.at(last);
	bool startOnConsonant = poem::vowelsRus[c];
	for (; last >= 0; last--) {
		unsigned char a = output.at(last);
		if (poem::vowelsRus[a]) break;
	}
	if (startOnConsonant) {
		if (last) {
			last--;
			c = output.at(last);
			if (c == (unsigned char)'ü'|| c == (unsigned char)'ú')
				last++;
			last += poem::vowelsRus[c];
		}
	}
	if (last < 0) return ""; // to eliminate dictionary errors
	string e;
	e.reserve(output.size() - last);

	for (int i = last; i < output.size(); i++) {
		unsigned char a = output.at(i);
		e += lettersToReplace[startOnConsonant][a];
	}
	return e;
}

void poem::updateRythm(int rhythmLenth, int rhythmRegularity) {
	if (rhythmRegularity < 2||rhythmLenth<2)
		return;
	if (rhythmLenth != syllableInfo.total || rhythmRegularity != syllableInfo.regularity) {
		syllableInfo.total = rhythmLenth;
		syllableInfo.regularity = rhythmRegularity;
		syllableInfo.firstAccent = (syllableInfo.total - 1) % syllableInfo.regularity;
		const wchar_t* nextParVocabularyFile = L"vocab_.txt";
		// Read word and slot index files.
		std::string str;
		size_t idx = 0;

		for (int i = 0; i < nodeType::rythmToIndex.size(); i++)
			for (int j = 0; j < nodeType::rythmToIndex[i].size(); j++)
				nodeType::rythmToIndex[i][j] = -1;
		//nodeType::rythmToIndex[0][0] = 0;
		poem::rythms.clear();
		vocabToIndex.clear();
		indexToSlots.clear();
		ifstream input(nextParVocabularyFile);
		while (getline(input, str)) {
			vocabToIndex[str] = idx++;
			nodeType n(str);
			indexToSlots.push_back(n);
		}
		/*std::ofstream out("output.txt");
		for (auto s : indexToSlots) {
			string ss;
			if (s.endIndex)
				ss = s.end();
			else
				ss = "";
			out << s.s <<" "<<ss<<"\r\n";
		}
		out.close();*/

		for (int i = 1; i < nodeType::rythmToIndex.size(); i++) { // every accentPos for zero or one syll is permitted
			nodeType::rythmToIndex[i][0] = nodeType::rythmToIndex[0][0];
			nodeType::rythmToIndex[i][1] = nodeType::rythmToIndex[0][1];
		}
		indexOf$ = vocabToIndex["$"];
		indexOfDot = vocabToIndex["."];
		indexOfSpace = vocabToIndex["~"];
		zeroMap[rythms[0][0]];
		text = "";
	}
}

poem::poem() {
	const wchar_t* nextParModel = L"lm_epoch0.dnn";

	// Load the model.
	nextParModelFunc= Function::Load(nextParModel, CNTK::DeviceDescriptor::CPUDevice());
	//insertSlotModelFunc = Function::Load(L"lm_insert.dnn", CNTK::DeviceDescriptor::CPUDevice());
	//replaceSlotModelFunc = Function::Load(L"lm_replace.dnn", CNTK::DeviceDescriptor::CPUDevice());
	//deleteSlotModelFunc = Function::Load(L"lm_delete.dnn", CNTK::DeviceDescriptor::CPUDevice());

	string s = "àå¸èîóûýþÿ";
	for (unsigned char a : s)
		vowelsRus[a] = true;
	for (int i = 0; i < 256; i++) {
		nodeType::lettersToReplace[0][i] = i;
		nodeType::lettersToReplace[1][i] = i;
	}
	string s1 = "áâãäæç÷";
	string s2 = "ïôêòøñù";
	for (unsigned int i = 0; i < s1.size(); i++) {
		unsigned char a = s1[i];
		nodeType::lettersToReplace[0][a] = s2[i];
		nodeType::lettersToReplace[1][a] = s2[i];
	}
	s1 = "ÿå¸èþ";
	s2 = "àýîûó";
	for (unsigned int i = 0; i < s1.size(); i++) {
		unsigned char a = s1[i];
		nodeType::lettersToReplace[0][a] = s2[i];
	}

	// Get input variables
	nextParInputVarForward = nextParModelFunc->Arguments()[0];
	nextParInputVarBackward = nextParModelFunc->Arguments()[1];
	//insertSlotInputVarForward = insertSlotModelFunc->Arguments()[0];
	//insertSlotInputVarBackward = insertSlotModelFunc->Arguments()[1];
	//replaceSlotInputVarForward = replaceSlotModelFunc->Arguments()[0];
	//deleteSlotInputVarForward = deleteSlotModelFunc->Arguments()[0];
	vocabSize = nextParInputVarForward.Shape().TotalSize();
	startSet[0] = 0;
	//syllableInfo = { 7,2,0 };
	updateRythm(8,2);
}

#include <regex>
#include "atlconv.h"
wstring string_to_wstring(string str) {
	wstring convertedString;
	int requiredSize = MultiByteToWideChar(CP_ACP, 0, str.c_str(), -1, 0, 0);
	if (requiredSize > 0) {
		vector<wchar_t> buffer(requiredSize);
		MultiByteToWideChar(CP_ACP, 0, str.c_str(), -1, &buffer[0], requiredSize);
		convertedString.assign(buffer.begin(), buffer.end() - 1);
	}
	return convertedString;
}

string wstring_to_string(wstring str) {
	char sAnsiStr[1024];
	int lenStr = WideCharToMultiByte(CP_ACP, 0, str.c_str(), -1, (char*)&sAnsiStr, 1024, NULL, NULL);
	string convertedString(sAnsiStr);
	return convertedString;
}

bool poem::updateText(string s) {
	if (s == text)
		return true;
	text = s;
	for (auto & c : s) // replace uppercase letters
		c = tolower(c);
	wstring ss = string_to_wstring(s.c_str());
	//ss = regex_replace(ss, wregex(L"$"), L"");// remove ends of line
	ss = regex_replace(ss, wregex(L"(\\w)- *(\\w)"), L"$1$2");// remove hyphens
	ss = regex_replace(ss, wregex(L"[0-9]+"), L"0");// merge digits
	ss = regex_replace(ss, wregex(L"[a-z]+"), L"a");// merge english
	ss = regex_replace(ss, wregex(L"[^\\w\\s`]+"), L"\.");// replace all punctuation by dots
	ss = regex_replace(ss, wregex(L"\\. *(?=\\.)"), L"");// merge dots
	ss = regex_replace(ss, wregex(L"^ *\\."), L"");// remove dot in the beginning
	ss = regex_replace(ss, wregex(L"(\\w)\\s*$"), L"$1.");// add dot to the end
	ss = regex_replace(ss, wregex(L"(\\w)\\."), L"$1 .");// add spaces after dots
	ss = regex_replace(ss, wregex(L"\\.(\\w)"), L". $1");// add spaces before dots
	s = wstring_to_string(ss.c_str());
	std::istringstream iss(s);
	std::vector<std::string> words((std::istream_iterator<std::string>(iss)),
		std::istream_iterator<std::string>());

	for (auto& w : words) {
		auto accPos = w.find('`');
		if (accPos == std::string::npos) {// add accents
			auto index = nodeType::accentizer.find(w);// if found in the dictionary - ok
			if (index == nodeType::accentizer.end()) {// else - calculate a number of vowels
				const std::string vowels = "àå¸èîóûýþÿ";  
				size_t n = 0, nvowels = 0, firstPos = 0;          
				while ((n = w.find_first_of(vowels, n)) != std::string::npos) {
					if (!nvowels)
						firstPos = n;
					nvowels++;  
					n++;        
				}
				if (nvowels > 1) {// if 2 or more vowels - accentize the first
					w.insert(firstPos, 1, '`');
				}
			}
			else {
				w = indexToSlots[index->second].s;
			}
		}
	}

	nodes.resize(words.size());
	for (int i = 0; i < words.size(); i++) {
		auto it= vocabToIndex.find(words[i]);
		if (it != vocabToIndex.end()) {
			nodes[i] = indexToSlots[it->second];
		}
		else {
			nodes[i] = nodeType(words[i]);
		}
	}

	wordsType::slotSeq.resize((2*nodes.size()-1)*vocabToIndex.size());

	// evaluate nextWord model
	for (int j = 0; j < 2; j++) {
		// create arrays of best words (priority queues) for every slot (not rythmed)
		std::vector<size_t> seqDataFwd;
		std::vector<size_t> seqDataBwd;
		seqDataFwd.push_back(indexOfDot);
		for (int i = 0; i < nodes.size() - 1-j; i++)
		{
			seqDataFwd.push_back(nodes[i].index);
		}
		for (int i = j; i < nodes.size(); i++)
		{
			seqDataBwd.push_back(nodes[i].index);
		}

		// Create input value using one-hot vector and input data map
		ValuePtr inputValForward = Value::CreateSequence<float>(vocabSize, seqDataFwd, CNTK::DeviceDescriptor::CPUDevice());
		ValuePtr inputValBackward = Value::CreateSequence<float>(vocabSize, seqDataBwd, CNTK::DeviceDescriptor::CPUDevice());
		std::unordered_map<Variable, ValuePtr> inputDataMap = { { nextParInputVarForward, inputValForward },{ nextParInputVarBackward, inputValBackward } };
		// The model has only one output.
		// If the model has more than one output, use modelFunc->Outputs to get the list of output variables.
		Variable outputVar = nextParModelFunc->Output();
		// Create output data map. Using null as Value to indicate using system allocated memory.
		// Alternatively, create a Value object and add it to the data map.
		std::unordered_map<Variable, ValuePtr> outputDataMap = { { outputVar, nullptr } };
		// Start evaluation on the device
		nextParModelFunc->Evaluate(inputDataMap, outputDataMap, CNTK::DeviceDescriptor::CPUDevice());
		// Get evaluate result as dense output
		ValuePtr outputVal = outputDataMap[outputVar];
		std::vector<std::vector<float>> outputData;
		outputVal->CopyVariableValueTo(outputVar, outputData);

		// output the result
		outputSampleSize = outputVar.Shape().TotalSize();
		if (outputData.size() != 1)
		{
			throw("Only one sequence of slots is expected as output.");
		}
		int numOfSlots = outputData[0].size() / outputSampleSize;
		for (int i = 0; i < numOfSlots; i++) {
			copy(outputData[0].begin()+i*outputSampleSize, outputData[0].begin() + (i+1) * outputSampleSize, wordsType::slotSeq.begin()+(i * 2 + j)*outputSampleSize);
		}
	}
	if (wordsType::slotSeq.size() % outputSampleSize != 0)
	{
		throw("The number of elements in the slot sequence is not a multiple of sample size");
	}

	// evaluate insert point model - useless - weak correllation
	/*vector <float> insertProbs(nodes.size()*2);
	fill(insertProbs.begin(), insertProbs.end(), 1);
	vector <float> deleteProbs(nodes.size());
	float totalReplace = 0;
	float totalDelete = 0;*/
	/*{
		// create arrays of insertion probabilities for every slot
		std::vector<size_t> seqDataFwd;
		std::vector<size_t> seqDataBwd;
		seqDataFwd.push_back(indexOfDot);
		for (int i = 0; i < nodes.size()-1; i++)
		{
			seqDataFwd.push_back(nodes[i].index);
		}
		for (int i = 0; i < nodes.size(); i++)
		{
			seqDataBwd.push_back(nodes[i].index);
		}

		// Create input value using one-hot vector and input data map
		ValuePtr inputValForward = Value::CreateSequence<float>(vocabSize, seqDataFwd, CNTK::DeviceDescriptor::CPUDevice());
		ValuePtr inputValBackward = Value::CreateSequence<float>(vocabSize, seqDataBwd, CNTK::DeviceDescriptor::CPUDevice());
		std::unordered_map<Variable, ValuePtr> inputDataMap = { { insertSlotInputVarForward, inputValForward },{ insertSlotInputVarBackward, inputValBackward } };
		// The model has only one output.
		// If the model has more than one output, use modelFunc->Outputs to get the list of output variables.
		Variable outputVar = insertSlotModelFunc->Output();
		// Create output data map. Using null as Value to indicate using system allocated memory.
		// Alternatively, create a Value object and add it to the data map.
		std::unordered_map<Variable, ValuePtr> outputDataMap = { { outputVar, nullptr } };
		// Start evaluation on the device
		insertSlotModelFunc->Evaluate(inputDataMap, outputDataMap, CNTK::DeviceDescriptor::CPUDevice());
		// Get evaluate result as dense output
		ValuePtr outputVal = outputDataMap[outputVar];
		std::vector<std::vector<float>> outputData;
		outputVal->CopyVariableValueTo(outputVar, outputData);

		// output the result
		int outputSampleSize = outputVar.Shape().TotalSize();
		if (outputData.size() != 1)
		{
			throw("Only one sequence of slots is expected as output.");
		}
		int numOfSlots = outputData[0].size() / outputSampleSize;
		float* a = &outputData[0][0];
		//float total = 0;
		for (int i = 0; i < numOfSlots; i++) {
			insertProbs[i*2] = exp(a[1]) / (exp(a[0]) + exp(a[1]));
			//total += insertProbs[i];
			a += 2;
		}
		if (outputData[0].size() % outputSampleSize != 0)
		{
			throw("The number of elements in the slot sequence is not a multiple of sample size");
		}
	}*/

	// evaluate replace point model
	/*{
		// create arrays of insertion probabilities for every slot
		std::vector<size_t> seqDataFwd;
		seqDataFwd.push_back(indexOfDot);
		for (int i = 0; i < nodes.size(); i++)
		{
			seqDataFwd.push_back(nodes[i].index);
		}

		// Create input value using one-hot vector and input data map
		ValuePtr inputValForward = Value::CreateSequence<float>(vocabSize, seqDataFwd, CNTK::DeviceDescriptor::CPUDevice());
		std::unordered_map<Variable, ValuePtr> inputDataMap = { { replaceSlotInputVarForward, inputValForward }};
		// The model has only one output.
		// If the model has more than one output, use modelFunc->Outputs to get the list of output variables.
		Variable outputVar = replaceSlotModelFunc->Output();
		// Create output data map. Using null as Value to indicate using system allocated memory.
		// Alternatively, create a Value object and add it to the data map.
		std::unordered_map<Variable, ValuePtr> outputDataMap = { { outputVar, nullptr } };
		// Start evaluation on the device
		replaceSlotModelFunc->Evaluate(inputDataMap, outputDataMap, CNTK::DeviceDescriptor::CPUDevice());
		// Get evaluate result as dense output
		ValuePtr outputVal = outputDataMap[outputVar];
		std::vector<std::vector<float>> outputData;
		outputVal->CopyVariableValueTo(outputVar, outputData);

		// output the result
		int outputSampleSize = outputVar.Shape().TotalSize();
		if (outputData.size() != 1)
		{
			throw("Only one sequence of slots is expected as output.");
		}
		int numOfSlots = outputData[0].size() / outputSampleSize;
		float* a = &outputData[0][2];
		for (int i = 0; i < numOfSlots-1; i++) {
			float p= exp(a[1]) / (exp(a[0]) + exp(a[1]));
			insertProbs[i * 2 + 1] = p;
			totalReplace += p;
			//insertProbs[i * 2 + 1] = a[1] / (a[0] + a[1]);
			a += 2;
		}
		if (outputData[0].size() % outputSampleSize != 0)
		{
			throw("The number of elements in the slot sequence is not a multiple of sample size");
		}
	}*/

	// evaluate delete point model
	/*{
		// create arrays of insertion probabilities for every slot
		std::vector<size_t> seqDataFwd;
		seqDataFwd.push_back(indexOfDot);
		for (int i = 0; i < nodes.size(); i++)
		{
			seqDataFwd.push_back(nodes[i].index);
		}

		// Create input value using one-hot vector and input data map
		ValuePtr inputValForward = Value::CreateSequence<float>(vocabSize, seqDataFwd, CNTK::DeviceDescriptor::CPUDevice());
		std::unordered_map<Variable, ValuePtr> inputDataMap = { { deleteSlotInputVarForward, inputValForward } };
		// The model has only one output.
		// If the model has more than one output, use modelFunc->Outputs to get the list of output variables.
		Variable outputVar = deleteSlotModelFunc->Output();
		// Create output data map. Using null as Value to indicate using system allocated memory.
		// Alternatively, create a Value object and add it to the data map.
		std::unordered_map<Variable, ValuePtr> outputDataMap = { { outputVar, nullptr } };
		// Start evaluation on the device
		deleteSlotModelFunc->Evaluate(inputDataMap, outputDataMap, CNTK::DeviceDescriptor::CPUDevice());
		// Get evaluate result as dense output
		ValuePtr outputVal = outputDataMap[outputVar];
		std::vector<std::vector<float>> outputData;
		outputVal->CopyVariableValueTo(outputVar, outputData);

		// output the result
		int outputSampleSize = outputVar.Shape().TotalSize();
		if (outputData.size() != 1)
		{
			throw("Only one sequence of slots is expected as output.");
		}
		int numOfSlots = outputData[0].size() / outputSampleSize;
		float* a = &outputData[0][2];
		for (int i = 0; i < numOfSlots - 1; i++) {
			float p= exp(a[1]) / (exp(a[0]) + exp(a[1]));
			deleteProbs[i] = p;
			totalDelete += p;
			a += 2;
		}
		if (outputData[0].size() % outputSampleSize != 0)
		{
			throw("The number of elements in the slot sequence is not a multiple of sample size");
		}
	}*/

	// create min heap
	pqWords[0] = priority_queue<wordsType>();
	pqWords[1] = priority_queue<wordsType>();
	size_t numOfSlotsInOutput = wordsType::slotSeq.size() / outputSampleSize;
	for (size_t slot = 0; slot < numOfSlotsInOutput; slot++) {
		int offset = slot * outputSampleSize;
		network1Output = &wordsType::slotSeq[slot * outputSampleSize];
		if (!(slot & 1) || nodes[slot / 2].index != indexOfDot) { // prohibits changing dots
			network1Output[indexOf$] = -100000000;
			network1Output[indexOfDot] = -100000000; // prohibits adding dots
			if (slot & 1) {
				network1Output[nodes[slot / 2].index] = -100000000;
			}
			float total = 0;
			for (int i = 0; i < outputSampleSize; i++) {
				float p = exp(network1Output[i]);
				total += p;
			}
			if (slot & 1) {
				network1Output[nodes[slot / 2].index] = log(total); // give the highest priority to initial words
				/*float p = exp(network1Output[0]);
				network1Output[0] = p / total*deleteProbs[slot/2]/totalDelete;
				wordsType w = { offset };
				pqWords[slot & 1].push(w);
				total = insertProbs[slot] / totalReplace / total;
				for (int i = 1; i < outputSampleSize; i++) {
					float p = exp(network1Output[i]);
					network1Output[i] = p * total;
					wordsType w = { offset + i };
					pqWords[slot & 1].push(w);
				}*/
			}
			/*else {
				total = insertProbs[slot] / total;
				for (int i = 0; i < outputSampleSize; i++) {
					float p = exp(network1Output[i]);
					network1Output[i] = p * total;
					wordsType w = { offset + i };
					pqWords[slot & 1].push(w);
				}
			}*/
			total = 1 / total;
			for (int i = 0; i < outputSampleSize; i++) {
				float p = exp(network1Output[i]);
				network1Output[i] = p * total;
				wordsType w = { offset + i };
				pqWords[slot & 1].push(w);
			}
		}
		else {
			network1Output[indexOfDot] = 1.f; 
			wordsType w = { offset + indexOfDot };
			pqWords[1].push(w);
		}
	}

	// add initial words to the top
	/*network1Output = &wordsType::slotSeq[0];
	float topWeight = wordsType::slotSeq[pqWords[1].top().number]+1;
	for (int slot = 0; slot < nodes.size()-1;slot++) {
		int ind = (2*slot+1) * outputSampleSize + nodes[slot].index;
		network1Output[ind] = topWeight;
		wordsType w = { ind };
		pqWords[1].push(w);
	}*/

	// Adding a word to slot with extraction of new sentences
	// It doesn't make sense to fill in all the above arrays based on the whole set of the words: it is impossible to
	// check based on syllToWord / endAndSyllToWord (rythm may exist, but the corresponding word is not added yet
	// The only option - try to update the arrays and check changes
	// in the very start - add words until no slots are empty (probably empty words to every slot)

	// create arrays
	syllAndEndToWords.clear();
	addedWords.clear();
	syllEndToWordsOrdered.clear();
	syllToWordsOrdered.clear();
	orderedWords.clear();
	wordToOrderedIndex.clear();
	endToSlotToSyll.clear();
	pairs.clear();
	accumulatedLengths.clear();
	for (auto&& f : fours)
		f.clear();

	syllAndEndToWords.resize(nodes.size() * 2-1);
	addedWords.resize(nodes.size()*2-1);

	syllEndToWordsOrdered.resize(nodes.size()*2-1);
	syllToWordsOrdered.resize(nodes.size()*2-1);

	nextDotPos.resize(nodes.size() * 2);
	int counter = nextDotPos.size();
	for (int i = nodes.size()-1; i >=0 ; i--) {
		if (nodes[i].index == indexOfDot) {
			counter = 2 * i + 1;
		}
		nextDotPos[2 * i+1] = counter;
		nextDotPos[2 * i] = counter;
	}

	prevDotPos.resize(nodes.size() * 2);
	counter = -1;
	for (int i = 0; i < nodes.size(); i++) {
		prevDotPos[2 * i] = counter;
		if (nodes[i].index == indexOfDot) {
			counter = 2 * i + 1;
		}
		prevDotPos[2 * i + 1] = counter;
	}
	nextDotPosPossible.resize(nodes.size() * 2);
	//prevDotPosPossible.resize(nodes.size() * 2);
	fill(nextDotPosPossible.begin(), nextDotPosPossible.end(), nextDotPosPossible.size());
	//fill(prevDotPosPossible.begin(), prevDotPosPossible.end(), -1);
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].index == indexOfDot) {
			int dp = 2 * i + 1;
			nextDotPosPossible[dp] = dp;
			//prevDotPosPossible[dp] = dp;
		}
	}
	spacePermitted.resize(nodes.size() * 2);
	fill(spacePermitted.begin(), spacePermitted.end(), false);
	zeroCurr.resize(nodes.size() * 2-1);
	zeroPrev.resize(nodes.size() * 2);
	zeroPrevPrev.resize(nodes.size() * 2 + 1);
	fill(zeroCurr.begin(), zeroCurr.end(), false);
	fill(zeroPrev.begin()+1, zeroPrev.end(), false);
	fill(zeroPrevPrev.begin()+2, zeroPrevPrev.end(), false);
	zeroPrev[0] = true;
	zeroPrevPrev[0] = true;
	zeroPrevPrev[1] = true;
	endSlotPermitted.resize(nodes.size() * 2);
	wordCounter = 0;

	for (int i = 0; i < nodes.size(); i++) {
		endSlotPermitted[2 * i] = rythmed(0, nodes[i].accentPos, nodes[i].syllableNumber);
		endSlotPermitted[2 * i + 1] = true;
	}
	nextRythmedPatternsPQ = priority_queue<patternType>();
	prevRythmedFours = stack<patternType>();
	return true;
}

string poem::nextWord() { // find rythmed pairs and then look through patterns
	//	when adding word into end slot: 
	//		when adding a pattern - calculate four end slots and keep them
	//		else for every pattern:
	//			for changed pair:
	//				compare n->endIndex against pairs[slot1][slot2][rythmPattern[syll1]][rythmPattern[syll2]]
	//				if not found: 
	//					check against syllAndEndToWords[oppositeSlot][rythmPattern[oppositeSlot]]
	//						if found: update pairs
	//			for opposite pair:
	//				check pairs[slot1][slot2][rythmPattern[syll1]][rythmPattern[syll2]]
	//			if both pairs exist - pattern is rythmed
	nextRythmedPatternsPQ = priority_queue<patternType>();
	prevRythmedFours = stack<patternType>();
	while (!(pqWords[0].empty()&& pqWords[1].empty())) {
		std::priority_queue< wordsType> *queue;
		int queueIndex;// = wordCounter < nodes.size() - 1 || wordCounter > 40 * (nodes.size() - 1);
		if (wordCounter < nodes.size() - 1)
			queueIndex = 1;
		else if (wordCounter < 40 * (nodes.size() - 1)) {
			queueIndex = 0;
			PQWeightRatio = wordsType::slotSeq[pqWords[0].top().number] / wordsType::slotSeq[pqWords[1].top().number];
		}
		else {
			//queueIndex = 1;
			//queueIndex = wordCounter & 1;
			queueIndex = wordsType::slotSeq[pqWords[0].top().number] < wordsType::slotSeq[pqWords[1].top().number] * PQWeightRatio;
		}
		wordCounter++;
		queue = &pqWords[queueIndex];
		addedWord = queue->top().number%outputSampleSize;
		slotOfNewWord = queue->top().number / outputSampleSize;
		nodeType* n;
		if (addedWord != indexOf$)
			n = &indexToSlots[addedWord];
		else
			n = &nodes[slotOfNewWord/2];
		rythmIndexOfAddedWord = n->rythmIndex;
		if (n->syllableNumber==0) {
			zeroCurr[slotOfNewWord] = true;
			zeroPrev[slotOfNewWord + 1] = true;
			zeroPrevPrev[slotOfNewWord + 2] = true;
			if (addedWord == indexOfSpace) { // update nextDotPos, prevDotPos
				spacePermitted[slotOfNewWord] = true;
				int ndp = nextDotPosPossible[slotOfNewWord + 1];// update nextDotPos
				if (ndp != nextDotPosPossible.size()) {
					nextDotPosPossible[slotOfNewWord] = ndp;
					for (int i = slotOfNewWord-1; i >= 0; i--) {
						if (spacePermitted[i])
							nextDotPosPossible[i] = ndp;
						else
							break;
					}
				}
				/*if (slotOfNewWord) {// update prevDotPos
					int pdp = prevDotPosPossible[slotOfNewWord - 1];
					if (pdp != -1) {
						prevDotPosPossible[slotOfNewWord] = pdp;
						for (int i = slotOfNewWord + 1; i < prevDotPosPossible.size(); i++) {
							if (spacePermitted[i])
								prevDotPosPossible[i] = pdp;
							else
								break;
						}
					}
				}*/
			}
		}
		endIndexOfAddedWord = n->endIndex;
		syllEndToWordsOrdered[slotOfNewWord][rythmIndexOfAddedWord][n->endIndex].push_back(orderedWords.size());
		syllToWordsOrdered[slotOfNewWord][rythmIndexOfAddedWord].push_back(orderedWords.size());
		wordToOrderedIndex[queue->top().number] = orderedWords.size();
		orderedWords.push_back(queue->top().number);
		queue->pop();
		bool updateFourNeeded = false;
		// check whether new end added
		if (n->endIndex!=nodeType::endStringToCode[""]&&endSlotPermitted[slotOfNewWord]) {
			// check whether another end already exists
			auto slotsForEnd = &endToSlotToSyll[n->endIndex];
			if (slotsForEnd->size()) {
				// check that end/slot/syll was not already checked
				auto syllsForSlot = slotsForEnd->find(slotOfNewWord);
				if (syllsForSlot != slotsForEnd->end()) {
					if (syllsForSlot->second.find(rythmIndexOfAddedWord) != syllsForSlot->second.end()) {
						continue;
					}
				}
				// look through opposite slot/sylls
				for (auto&& opp : *slotsForEnd) {
					if (abs(opp.first - slotOfNewWord) < 2)
						continue;
					for (auto syll : opp.second) {
						vector<vector<int>> p;
						if (slotOfNewWord<opp.first)
							p={ { slotOfNewWord, opp.first },{ rythmIndexOfAddedWord, syll } };
						else
							p={ { opp.first, slotOfNewWord },{ syll, rythmIndexOfAddedWord } };
						// add to fours
						// look through other pairs
						for (auto&& p2sl0 : pairs) {
							bool firstIsLess0 = p[0][0] < p2sl0.first && p[0][1] > p2sl0.first;
							bool secondIsLess0 = !firstIsLess0&&p[0][0] > p2sl0.first;
							if (firstIsLess0 || secondIsLess0) {
								for (auto&& p2sl1 : p2sl0.second) {
									bool firstIsLess = firstIsLess0&&p[0][1] < p2sl1.first;
									bool secondIsLess = secondIsLess0&&p[0][1] > p2sl1.first && p[0][0] < p2sl1.first;
									if (firstIsLess || secondIsLess) {
										int ind1 = firstIsLess;
										int ind2 = !firstIsLess;
										// add four
										if (firstIsLess) {
											slots = { -1,p[0][0],p2sl0.first,p[0][1],p2sl1.first };
										}
										else {
											slots = { -1,p2sl0.first,p[0][0],p2sl1.first,p[0][1] };
										}
										int endSlot = slots[4];
										if (!(endSlot & 1) && endSlot + 2<syllAndEndToWords.size() && nodes[endSlot / 2].syllableNumber) {
											continue;// break;
										}
										// create/extend accLen tables
										for (int lineIndex = 0; lineIndex < 4; lineIndex++) {
											int startSlot = slots[lineIndex] + 1;
											int endSlot = slots[lineIndex + 1];
											vector<array<map<int,int>, 2>> *al;
											int alSize;
											if (startSlot&1 && !rythmed(0, nodes[startSlot / 2].accentPos, nodes[startSlot / 2].syllableNumber)) {
												al=NULL;
												alSize = 0;
											}
											else {
												al = &accumulatedLengths[startSlot];
												alSize = al->size();
											}
											if (endSlot - startSlot >= alSize) {
												if(al) 
													updateAccumLengthsForNewEnd(al, startSlot, endSlot);
												// create/extend accLen tables for other scenarios
												int firstDotPos = nextDotPos[startSlot];
												bool dot = firstDotPos < endSlot;
												if (dot) {// update table for the second scenario
													startSlot = firstDotPos + 1;
													al = &accumulatedLengths[startSlot];
													if (endSlot - startSlot >= al->size()) {
														updateAccumLengthsForNewEnd(al, startSlot, endSlot);
														int lastDotPos = prevDotPos[endSlot];
														bool dot2 = lastDotPos > firstDotPos + 1;
														if (dot2) { // update table for the third scenario
															startSlot = lastDotPos + 1;
															al = &accumulatedLengths[startSlot];
															if (endSlot - startSlot >= al->size()) {
																updateAccumLengthsForNewEnd(al, startSlot, endSlot);
															}
														}
													}
												}
											}
										}
										for (auto&& p2syll0 : p2sl1.second) {
											for (auto&& p2syll1 : p2syll0.second) {
												if (firstIsLess) {
													sylls = { 0,p[1][0],p2syll0.first,p[1][1],p2syll1.first };
												}
												else {
													sylls = { 0,p2syll0.first,p[1][0],p2syll1.first,p[1][1] };
												}
												pairs[p[0][0]][p[0][1]][p[1][0]][p[1][1]].insert(n->endIndex);
												updateFourForNewEndRecoursively(0, true);
												updateFourNeeded = true;
											nextPair:;
											}
										}
									}
								}
							}
						}
						pairs[p[0][0]][p[0][1]][p[1][0]][p[1][1]].insert(n->endIndex);
					}
				}
			}
			(*slotsForEnd)[slotOfNewWord].insert(rythmIndexOfAddedWord);
		}
		// check whether new syll added
		if (syllAndEndToWords[slotOfNewWord].find(rythmIndexOfAddedWord) == syllAndEndToWords[slotOfNewWord].end()) {
			syllAndEndToWords[slotOfNewWord][rythmIndexOfAddedWord][n->endIndex] = { addedWord }; // add to syllAndEndToWord
			// new syll added - update all accLenths
			for (auto&& slot0 : accumulatedLengths) {
				// update accLen if the subpattern changed
				int startSlotOfALVector = slot0.first;
				auto al = &accumulatedLengths[slot0.first]; // vector
				if (slotOfNewWord >= startSlotOfALVector && slotOfNewWord - startSlotOfALVector < al->size()) {
					int locSlot = slotOfNewWord - startSlotOfALVector;
					map<int,int>*  prevSet;
					bool lastAlti = startSlotOfALVector;
					int firstAlti = startSlotOfALVector & 1;
					for (int alti = firstAlti; alti <= lastAlti; alti++) { // accLenTableIndex
						if (!locSlot) {
							prevSet = &startSet;
						}
						else {
							prevSet = &(*al)[locSlot - 1][alti];
						}
						if (!locSlot) {
							if (alti) {
								if (slotOfNewWord & 1) {
									continue;
								}
								else {
									if (rythmIndexOfAddedWord) {
										continue;
									}
								}
							}
						}
						updateAccumLengthsForNewSyllRecoursively(al, slotOfNewWord, locSlot, prevSet, alti);
					}
				}
			}
			// look through fours recoursively
			updateFourNeeded = true;
		}
		else
			syllAndEndToWords[slotOfNewWord][rythmIndexOfAddedWord][n->endIndex].insert(addedWord);// add to syllAndEndToWord
		if (updateFourNeeded) {
			updateFourForNewSyllRecoursively(0, 0, true);
		}
		if (nextRythmedPatternsPQ.size()) {
			return bestSentenceOfBestFour();
		}
	}
	return "";
}

bool poem::accLenIsPermitted(vector<array<map<int, int>, 2>> *al, int syll, int prevAccLen, int locSlot, int accLenTableIndex, int limitedConsistency) {
	int len = rythms[syll][1];
	int newAccLen = prevAccLen + len;
	if (newAccLen > syllableInfo.total)
		return false;
	bool permitted = rythmed(prevAccLen, rythms[syll][0], len);
	if (permitted) {
		// how to prevent insertion of zero accLen after dot in first scenario ?
		// it will be prevented iff no zero accLen in the slot of dot
		auto s = &(*al)[locSlot][accLenTableIndex];
		auto l = s->find(newAccLen);
		if (l == s->end() || (!limitedConsistency && l->second)) {
			(*s)[newAccLen]= limitedConsistency;
			if (locSlot < al->size() - 1)
				return true;
		}
	}
	return false;
}

void poem::updateAccumLengthsForNewSyllRecoursively(vector<array<map<int, int>, 2>> *al, int slot, int locSlot, const map<int, int>* prevSet, int alti) {
	if (slot & 1) {
		int originalSyll = nodes[slot / 2].rythmIndex;
		// first - create set for original syll and mark them by 0
		for (auto prev : *prevSet) {
			if (accLenIsPermitted(al, originalSyll, prev.first, locSlot, alti,0)) {// recourse;
				updateAccumLengthsForNewSyllRecoursively(al, slot + 1, locSlot + 1, &(*al)[locSlot][alti], alti);
			}
		}
		// second - if zeroPrev - check sylls other than original and mark them by 1
		if (locSlot || !alti) {
			if (zeroPrev[slot]) {
				for (auto&& syll : syllAndEndToWords[slot]) {
					if (syll.first != originalSyll) {// add prevPrev+syll
						for (auto prev : *prevSet) {
							if (!prev.second) {// if limitedConsistency - not all prev should be used
								if (accLenIsPermitted(al, syll.first, prev.first, locSlot, alti, 1)) {// recourse;
									updateAccumLengthsForNewSyllRecoursively(al, slot + 1, locSlot + 1, &(*al)[locSlot][alti], alti);
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		// first - create set for zero syll and mark them by 0
		if (zeroCurr[slot]) {
			for (auto prev : *prevSet) {
				if (accLenIsPermitted(al, rythms[0][0], prev.first, locSlot, alti,0)) {// recourse;
					updateAccumLengthsForNewSyllRecoursively(al, slot + 1, locSlot + 1, &(*al)[locSlot][alti], alti);
				}
			}
		}
		// second - check sylls other than zero and mark them by 1
		if (locSlot || !alti) {
			for (auto&& syll : syllAndEndToWords[slot]) {
				if (syll.first != rythms[0][0]) {
					for (auto prev : *prevSet) {
						if (!prev.second) {// if limitedConsistency - not all prev should be used
							if (accLenIsPermitted(al, syll.first, prev.first, locSlot, alti, 1)) {// recourse;
								updateAccumLengthsForNewSyllRecoursively(al, slot + 1, locSlot + 1, &(*al)[locSlot][alti], alti);
							}
						}
					}
				}
			}
		}
	}
}

void poem::fillNextAccLen(map<int, int> *prevSet, int syll, map<int, int> *res, int limitedConsistency) {
	int len = rythms[syll][1];
	int ap = rythms[syll][0];
	for (auto prev : *prevSet) {
		if (!(limitedConsistency && prev.second)) {// if limitedConsistency - not all prev should be used
			// check that the accLen doesn't exceed maximum
			int newAccLen = prev.first + len;
			if (newAccLen > syllableInfo.total)
				break; // acclengths are sorted 
			bool permitted = rythmed(prev.first, ap, len);
			if (permitted) {
				if (res->find(newAccLen) == res->end()) {
					(*res)[newAccLen] = limitedConsistency;
				}
			}
		}
	}
}

void poem::updateAccumLengthsForNewEnd(vector<array<map<int, int>, 2>> *al, int startSlotOfALVector, int endSlot) {
	assert(al->size() < endSlot - startSlotOfALVector + 1);
	int startSlot = startSlotOfALVector + al->size();
	al->resize(endSlot - startSlotOfALVector + 1);
	map<int, int>* prevSet; // prev, prevPrev, prevPrevPrev
	int startAlti, endAlti;
	if (startSlotOfALVector) 
		if (startSlotOfALVector & 1) {
			startAlti = 1;
			endAlti = 1;
		}
		else {
			startAlti = 0;
			endAlti = 1;
		}
	else {
		startAlti = 0;
		endAlti = 0;
	}
	for (int alti=startAlti; alti <= endAlti; alti++) {
		int locSlot = startSlot - startSlotOfALVector;
		if (!locSlot) {
			prevSet = &startSet;
		}
		else {
			prevSet = &(*al)[locSlot - 1][alti];
		}
		for (int slot = startSlot; slot <= endSlot; slot++) {
			auto alSlot = &(*al)[locSlot][alti];
			if (slot & 1) {
				// first - create set for original syll and mark them by 0
				int originalSyll = nodes[slot / 2].rythmIndex;
				fillNextAccLen(prevSet, originalSyll, alSlot,0);// add prev+syll
				// second - if zeroPrev - check sylls other than original and mark them by 1
				if (zeroPrev[slot]) {
					if (locSlot || !alti) {
						for (auto&& syll : syllAndEndToWords[slot]) {
							if (syll.first != originalSyll) {// if differs from original - start from prevPrev and fill the next slot
								fillNextAccLen(prevSet, syll.first, alSlot,1);// add prevPrev+syll
							}
						}
					}
				}
			}
			else {
				// first - create set for zero syll and mark them by 0
				if (zeroCurr[slot]) {
					fillNextAccLen(prevSet, rythms[0][0], alSlot,0);// add prev+syll
				}
				// second - check sylls other than zero and mark them by 1
				if (locSlot || !alti) {
					for (auto&& syll : syllAndEndToWords[slot]) {
						if (syll.first != rythms[0][0]) {
							fillNextAccLen(prevSet, syll.first, alSlot,1);// add prev+syll
						}
					}
				}
			}
			if (!alSlot->size()) {
				break;
			}
			prevSet = alSlot;
			locSlot++;
		}
	}
}

void poem::updateFourForNewSyllRecoursively(const int lineIndex, const int slot, const bool prevLineHaveNoWordBeforeDot) {
	/*counter++;
	if (counter == 8215)
		float ggg = 0;*/
	int altiForScen[3] = { accLenTableIndex(lineIndex) ,0,0 };
	int maxScen;
	if (altiForScen[0] && (slot & 1) && nodes[slot / 2].syllableNumber)
		maxScen = 1;
	else
		maxScen = 3;
	// look through fours and check whether they are permitted
	auto nextSlots = &fours[lineIndex][slot];
	for (auto slot1 = nextSlots->begin(); slot1 != nextSlots->end();slot1++) {
		int firstDotPos = nextDotPosPossible[slot];
		bool dot = firstDotPos < slot1->first;
		//bool dot2 = lastDotPos > firstDotPos + 1;
		int startSlot = slot;
		for (int scenario = 0; scenario < maxScen; scenario++) {
			if (dot) {
				// change the first slot for scenario 1 and 2
				if (scenario == 1) {
					startSlot = firstDotPos + 1;
				}
				else if (scenario == 2) {
					// check that the last dot is accessible from the first dot
					int lastDotPos = firstDotPos;// prevDotPosPossible[slot1->first];
					int pdp = prevDotPos[slot1->first];
					while (lastDotPos<pdp) {
						lastDotPos = nextDotPosPossible[lastDotPos + 1];
						if (lastDotPos == pdp) {
							break;
						}
					}
					if (lastDotPos != pdp)
						break;
					startSlot = lastDotPos + 1;
					if (startSlot <= firstDotPos + 1) {
						break;
					}
					//float fff = 0;
					/*int slot = firstDotPos;
					do {
						slot = nextDotPosPossible[slot+1];
					} while (slot < lastDotPos);
					if (slot != lastDotPos)
						break;*/
				}
				else if (lineIndex) {// dont check the first scenario if dot and not first line
					continue;// no need to check next slots (sorted)
				}
			}
			else {
				if (scenario || !prevLineHaveNoWordBeforeDot) {// prohibited (scenario for the second line only)
					break;// next slots should be checked (may include dot)
				}
			}
			auto al = &accumulatedLengths[startSlot]; // vector
			assert(al->size());
			int endSlot = slot1->first - startSlot;
			int alti = altiForScen[scenario];
			for (auto syll1 = slot1->second.begin(); syll1 != slot1->second.end();syll1++) {
				bool onlyNonOriginalEnd=false;
				if (lineIndex > 1) {// check that the pair exists for lines 2 and 3
					int i = (lineIndex & 1) + 1;
					auto sl1 = pairs.find(slots[i]);
					if (sl1 != pairs.end()) {
						auto sl2 = sl1->second.find(slot1->first);
						if (sl2 != sl1->second.end()) {
							auto syl1 = sl2->second.find(sylls[i]);
							if (syl1 != sl2->second.end()) {
								auto syl2 = syl1->second.find(syll1->first);
								if (syl2 == syl1->second.end()) {
									continue;
								}
								else {
									// check that the opposite is permitted
									onlyNonOriginalEnd = slots[i] & 1 && syl2->second.find(nodes[slots[i] / 2].endIndex) == syl2->second.end();
									if (onlyNonOriginalEnd&&!alternativeEndsPermitted[i-1]) {
										continue;
									}
									// check that the current is permitted
									onlyNonOriginalEnd = slot1->first & 1 && syl2->second.find(nodes[slot1->first / 2].endIndex) == syl2->second.end();
								}
							}
							else {
								goto nextSlot;
							}
						}
						else {
							goto nextSlot;
						}
					}
					else {
						return;
					}
				}
				int len = rythms[syll1->first][1];
				bool recourse = false;
				bool altEndPermitted = false;
				// if within slots and rythm is false - check whether the newSlot is rythmed 
				for (auto endScenario = syll1->second.begin(); endScenario != syll1->second.end();endScenario++) {// look through endScenarios
					if (!endScenario->first&&onlyNonOriginalEnd)
						continue;
					int nextSlot = endScenario->second[scenario];
					bool permitted=false;// = nextSlot >= 0;
					if (nextSlot!=-2) {
						// check whether the subpattern permitted
						//if (slotOfNewWord >= startSlot && slotOfNewWord <= slot1->first) {
							int delta, deltaLen;
							if (!(slot1->first & 1)) {
								nodeType* n = &nodes[(slot1->first - 1) / 2];
								//int accLenPrevPrev = syllableInfo.total - len - n->syllableNumber;
								delta = 1;
								deltaLen = len + n->syllableNumber;
							}
							else {
								deltaLen = len;
								if (endScenario->first)
									delta = 1;
								else {
									bool lastSyllChanged = syll1->first != nodes[slot1->first / 2].rythmIndex;
									delta = lastSyllChanged;
								}
							}
							if (endSlot - delta>0) {
								auto s = &(*al)[endSlot - 1 - delta][alti]; // set
								permitted = ((*s).find(syllableInfo.total - deltaLen) != s->end());
							}
							else {
								if (endSlot == delta) {
									if (endSlot) {
										auto s = &(*al)[endSlot - 1][alti]; // set
										permitted = deltaLen == syllableInfo.total && ((*s).find(syllableInfo.total - len) != s->end());
									}
									else {
										permitted = len == syllableInfo.total;
									}
								}
								else {
									permitted = len == syllableInfo.total && !alti;
								}
							}
							/*if (permitted) {
								endScenario->second[scenario] = startSlot;
							}*/
						//}
					}
					recourse |= permitted;
					altEndPermitted = !endScenario->first && permitted;
				}
				if (recourse) { // recourse and check further
					slots[lineIndex + 1] = slot1->first;
					sylls[lineIndex + 1] = syll1->first;
					currentScenarios[lineIndex] = startSlot;
					if (lineIndex < 3) {
						if(lineIndex<2)
							alternativeEndsPermitted[lineIndex & 1] = altEndPermitted;
						int firstDotPos = nextDotPos[slot];
						bool dot = firstDotPos < slot1->first;
						int lastDotPos = prevDotPos[slot1->first];
						bool dot2 = lastDotPos > firstDotPos + 1;
						updateFourForNewSyllRecoursively(lineIndex + 1, slot1->first + 1, (scenario || !dot)&&(scenario!=1||!dot2));
					}
					else {
						// check that there are zero syll after the last end
						for (int slot = slots[4] + 1; slot < syllAndEndToWords.size(); slot++) {
							if (!zeroCurr[slot])
								goto nextEndScenario;
						}
						// add patterns based on slots and sylls
						scenarios.push_back(currentScenarios);
						findRythmedPatternsBasedOnFour();
					nextEndScenario:;
					}
				}
			}
		}
	nextSlot:;
	}
}

int poem::accLenTableIndex(int lineIndex) {
	if (lineIndex) {
		if (slots[lineIndex] & 1) { // prev last syll is original or changed
			return sylls[lineIndex] != nodes[slots[lineIndex] / 2].rythmIndex;
		}
		else { // prev last syll is zero or not
			return 1;// !rythms[sylls[lineIndex]][1];
		}
	}
	else {
		return 0;
	}
}

bool poem::updateFourForNewEndRecoursively(int lineIndex, bool prevLineHaveNoWordBeforeDot) {
	// all lines should be added (no lines added in a course of updateFourForNewSyll
	// dont check accLenghts here.
	// it is not enough to check accLen table on the end of line - it could be unaccessible for the last syll
	// under certain conditions - end syll from the prev line is important. 
	// Possible scenarios: lastSyllChanged/Origin and lastSyllZero/NotZero. It limits the first syll:
	// if(slots[lineIndex]&1) 
	//    if(lastSyllChanged)
	//       only zero syll permitted in slots[lineIndex]+1
	// else if(lastSyllNotZero)
	//       only original syll permitted in slots[lineIndex]+1
	// But it is impossible to use with accLen tables (two separate accLen table needed for different scenarios)
	int endSlot = slots[lineIndex + 1];
	bool success0, // start from beginning
		success1, // start from first dot
		success2; // start from last dot
	int firstDotPos = nextDotPos[slots[lineIndex] + 1];
	bool dot = firstDotPos < endSlot;
	if (!dot && !prevLineHaveNoWordBeforeDot) {// prohibited
		return false;
	}
	int lastDotPos = prevDotPos[endSlot];
	bool dot2 = lastDotPos > firstDotPos + 1;
	bool fourPossible; // is the pattern possible in the future after new sylls added
	int delta, deltaLen;
	int len = rythms[sylls[lineIndex + 1]][1];
	// determine accLen table to use (common or limited)
	int accLenTI = accLenTableIndex(lineIndex); // 0 -common, 1- limited
	int endSlotLoc = endSlot - slots[lineIndex] - 1;
	if (!(endSlot & 1)) {
		if (!endSlotLoc) {
			fourPossible = len==syllableInfo.total&&!accLenTI;
			deltaLen = len;
		}
		else {
			nodeType* n = &nodes[(endSlot - 1) / 2];
			int accLenPrevPrev = syllableInfo.total - len - n->syllableNumber;
			fourPossible = rythmed(accLenPrevPrev, n->accentPos, n->syllableNumber);
			deltaLen = len + n->syllableNumber;
			if (endSlotLoc == 1) {
				fourPossible = !accLenPrevPrev;
			}
		}
		delta = 1;
	}
	else {
		deltaLen = len;
		bool lastSyllChanged = len != nodes[endSlot / 2].syllableNumber;
		delta = lastSyllChanged;
		if (endSlotLoc == 0) {
			fourPossible = !lastSyllChanged&&len == syllableInfo.total;
		}
		else {
			if (endSlotLoc == 1 && lastSyllChanged) {
				fourPossible = len == syllableInfo.total;
			}
			else {
				fourPossible = true;
			}
		}
	}
	if (!fourPossible)
		return false;
	auto m = &fours[lineIndex][slots[lineIndex] + 1][endSlot];
	if (m->find(sylls[lineIndex + 1]) == m->end()) { // new subpattern
													 // check that it is permitted
		int startSlot = slots[lineIndex] + 1;
		bool onlyScen0Possible = accLenTI && (startSlot & 1) && nodes[startSlot / 2].syllableNumber;
		// determine possible endScenarios (not only for the added word - because the procedure is recoursive)
		int maxEndScenario = endSlot & 1 && 
			syllAndEndToWords[endSlot][sylls[lineIndex + 1]].size()>1; // check that the alternative end is already added
		success0 = !dot || !lineIndex;
		success1 = !onlyScen0Possible && dot;
		success2 = !onlyScen0Possible && dot2;
		for (int endScenario = 0; endScenario <= maxEndScenario; endScenario++) {
			(*m)[sylls[lineIndex + 1]][endScenario] = { -2 + success0,-2 + success1,-2 + success2 };
		}
	}
	else {
		// determine possible endScenarios (not only for the added word - because the procedure is recoursive)
		auto it = (*m)[sylls[lineIndex + 1]].begin();
		auto dotPositions = it->second;
		success0 = dotPositions[0] >= -1;
		success1 = dotPositions[1] >= -1;
		success2 = dotPositions[2] >= -1;
	}
	if (lineIndex < 3) {
		// check the second dot - doesn't make sense - the scheme is not appropriate 
		// two scenarios are not enough in the case of two dots - third is needed - zeros before the second dot
		if (success0) {
			success0 = updateFourForNewEndRecoursively(lineIndex + 1, !dot);
		}
		if (success1) {
			success1 = updateFourForNewEndRecoursively(lineIndex + 1, !dot2);
		}
		if (success2) {
			success2 = updateFourForNewEndRecoursively(lineIndex + 1, true);
		}
	}
	return success0 || success1 || success2; // reports whether the full four succeeded
}

void poem::findRythmedPatternsBasedOnFour() {
	// look through all combinations of dims
	// how to compare different fours ? It is impossible to merge different fours in one graph (end slots should be paired).
	// The only way - using queue of queues. Weights of best ends should be added
	patternType pattern;
	for (int i = 1; i < 3; i++) {// check that the words are different
		auto&& ends1 = pairs[slots[i]][slots[i+2]][sylls[i]][sylls[i + 2]];
		for (auto&& e1 : ends1) {
			auto&& words1 = syllAndEndToWords[slots[i]][sylls[i]][e1];
			auto&& words2 = syllAndEndToWords[slots[i + 2]][sylls[i + 2]][e1];
			if (words1.size() == 1 && words2.size() == 1 && *words1.begin() == *words2.begin())
				continue;
			// add end to pattern
			pattern.ends[i-1].insert(e1);
		}
	}
	if (pattern.ends[0].size() && pattern.ends[1].size()) {
		//pattern.ends = { pairs[slots[1]][slots[3]][sylls[1]][sylls[3]],pairs[slots[2]][slots[4]][sylls[2]][sylls[4]] };
		//assert(pattern.ends[0].size() && pattern.ends[1].size());
		for (auto currentDimScenarios : scenarios) {
			for (int lineIndex = 0; lineIndex < 4; lineIndex++) {
				// remove doubling - exclude zero syll from accLen on dot positions
				int dotPos; // dont add zero accLen to avoid doubling
				int slot2 = slots[lineIndex] + 1;
				int firstDotPos = nextDotPos[slot2];
				int currentDimScenario = currentDimScenarios[lineIndex];
				if (currentDimScenario != slot2) {
					if (currentDimScenario == firstDotPos + 1) { // second scenario - remove zero accLen from lastDotPos
						int lastDotPos = prevDotPos[slots[lineIndex + 1]];
						if (lastDotPos > firstDotPos + 1) {
							dotPos = lastDotPos - slot2;
						}
						else {
							dotPos = -1;
						}
					}
					else {
						dotPos = -1;
					}
				}
				else { // first scenario - remove zero accLen from firstDotPos
					dotPos = firstDotPos - slot2;
				}
				vector<set<int>> accLen;
				accLen.resize(slots[lineIndex + 1] - slots[lineIndex]);
				int slot = 0;
				for (; slot < currentDimScenario - slot2; slot++) {// copy from prev
					if (slot) {
						accLen[slot] = accLen[slot - 1];
					}
					else {
						accLen[slot].insert(0);
					}
				}
				auto al = &accumulatedLengths[currentDimScenario]; // vector
				int maxLen = syllableInfo.total - rythms[sylls[lineIndex + 1]][1];
				int delta = slot;
				int alti;
				if (!delta)
					alti = accLenTableIndex(lineIndex);
				else alti = 0;
				do {
					if (slot == slots[lineIndex + 1] - slot2) {
						accLen[slot].insert(syllableInfo.total);
						break;
					}
					else {
						if (slot == slots[lineIndex + 1] - 1 - slot2) {
							accLen[slot].insert(maxLen);
						}
						else {
							for (auto l : (*al)[slot - delta][alti]) { // up to maxLen
								if (l.first > maxLen)
									break;
								// remove zero accLen in dot positions when needed
								if (slot != dotPos || l.first>0)
									accLen[slot].insert(l.first);
							}
						}
					}
					slot++;
				} while (true);
				vector<set<int>::iterator> patternIterator;
				patternIterator.resize(accLen.size());
				// init currentPattern
				for (int i = 0; i <= delta; i++) {
					patternIterator[i] = accLen[i].begin();
				}
				slot = delta;
				int prevAccLen = 0;
				int syll;
				vector<int> rythmPatternLoc(accLen.size());
				do {
					do {
						while (patternIterator[slot] == accLen[slot].end()) {
							if (slot) {
								slot--;
								patternIterator[slot]++;
								if (slot)
									prevAccLen = *patternIterator[slot - 1];
								else
									prevAccLen = 0;
							}
							else {
								goto endOfCycle;
							}
						}
						int syllLen = *patternIterator[slot] - prevAccLen;
						if (syllLen < nodeType::rythmToIndex[0].size()) {
							int accentPos = (syllableInfo.regularity-(prevAccLen % syllableInfo.total + syllableInfo.regularity - syllableInfo.firstAccent) % syllableInfo.regularity) % syllableInfo.regularity;
							syll = nodeType::rythmToIndex[accentPos][syllLen];
							if (syllAndEndToWords[slot + slot2].find(syll) == syllAndEndToWords[slot + slot2].end()) {
								patternIterator[slot]++;
							}
							else {
								if ((slot2 + slot) & 1) {
									int prevSyll;
									if (slot) {
										prevSyll = rythmPatternLoc[slot - 1];
									}
									else { // check the last word in the prev line
										prevSyll = sylls[lineIndex];
									}
									if (prevSyll) {// if previous is not zero - only original rythm is permitted
										if (syllLen != nodes[(slot + slot2) / 2].syllableNumber) {
											patternIterator[slot]++;
											continue;
										}
									}
								}
								else {
									if (slot2 + slot) {
										int prevSyll;
										if (slot) {
											prevSyll = rythmPatternLoc[slot - 1];
										}
										else { // check the last word in the prev line
											prevSyll = sylls[lineIndex];
										}
										if (prevSyll != nodes[(slot + slot2 - 1) / 2].rythmIndex) {// if previous is not original rythm - only zero is permitted
											if (syllLen) {
												patternIterator[slot] = accLen[slot].end(); // break because of the set is sorted
												continue;
											}
										}
									}
								}
								rythmPatternLoc[slot] = syll;
								slot++;
								if (slot < patternIterator.size()) {
									prevAccLen = *patternIterator[slot - 1];
									patternIterator[slot] = accLen[slot].lower_bound(prevAccLen);
								}
								else
									break;
							}
						}
						else {
							patternIterator[slot] = accLen[slot].end(); // break because of the set is sorted
						}
					} while (true);
				endOfCycle:;
					if (patternIterator[0] != accLen[0].end()) { // add pattern to queue
						vector<int> weights(accLen.size() - 1); // don't check last word (weight is determined by end)
						for (int pos = 0; pos < weights.size(); pos++) { // calculate preliminary weight based on sylls
							int rp = rythmPatternLoc[pos];
							weights[pos] = syllToWordsOrdered[pos + slot2][rp][0];
						}
						sort(weights.begin(), weights.end(), greater<int>());
						// check that the pattern was not already added
						pattern.nextRythmedPatternsPQ[lineIndex].push({ rythmPatternLoc,weights });
						patternIterator.back()++;
						slot = patternIterator.size() - 1;
						if (slot) {
							prevAccLen = *patternIterator[slot - 1];
						}
						else {
							prevAccLen = 0;
						}
					}
					else {
						//assert(pattern.nextRythmedPatternsPQ[lineIndex].size()>0);
						// prohobiting zero accLength on dot positions may rarely prohibit pattern when pattern is possible only 
						// if accLen == 0 on dot some point. Such situation can not be found. 
						if (!pattern.nextRythmedPatternsPQ[lineIndex].size()) {
							scenarios.clear();
							return;
						}
						break;
					}
				} while (true);
			}
		}
		// calculate pattern weight and push to queue
		for (int i = 0; i < 4; i++) {
			pattern.weight.insert(pattern.weight.end(), pattern.nextRythmedPatternsPQ[i].top().weight.begin(), pattern.nextRythmedPatternsPQ[i].top().weight.end()); // append to pattern weights
			int end = *pattern.ends[i & 1].begin();
			int endWeight = syllEndToWordsOrdered[slots[i + 1]][sylls[i + 1]][end][0];
			pattern.weight.push_back(endWeight); // append weight of the best end
		}
		// add last weights (after slots[4])
		for (int pos = slots[4] + 1; pos < syllToWordsOrdered.size(); pos++) { // calculate preliminary weight based on sylls
			pattern.weight.push_back(syllToWordsOrdered[pos][rythms[0][0]][0]);
		}
		sort(pattern.weight.begin(), pattern.weight.end(), greater<int>());
		nextRythmedPatternsPQ.push(pattern);
	}
	scenarios.clear();
}

string poem::bestSentenseForPattern() {
	if (!nextRythmedPatternsPQ.size())
		return "";
	rythmPattern.clear();
	for (int lineIndex = 0; lineIndex < 4; lineIndex++) {
		rythmPattern.insert(rythmPattern.end(), currentPattern.nextRythmedPatternsPQ[lineIndex].top().sylls.begin(), currentPattern.nextRythmedPatternsPQ[lineIndex].top().sylls.end());
		slots[lineIndex + 1] = rythmPattern.size() - 1;
	}
	/*for (int i = rythmPattern.size(); i < syllAndEndToWords.size(); i++)
		rythmPattern.push_back(0);*/
	return bestSentence();
}

string poem::bestSentenceOfBestFour() {
	if (!nextRythmedPatternsPQ.size())
		return "";
	currentPattern = nextRythmedPatternsPQ.top();
	ends = { *currentPattern.ends[0].begin(),*currentPattern.ends[1].begin() };
	return bestSentenseForPattern();
}

set<string> poem::possibleEnds(int index) {
	set<string> res;
	if (nextRythmedPatternsPQ.size()) {
		for (auto e : nextRythmedPatternsPQ.top().ends[index]) {
			res.insert(nodeType::endCodeToString[e]);
		}
	}
	return res;
}

array<string, 2> poem::currentEnds() {
	return { nodeType::endCodeToString[ends[0]],nodeType::endCodeToString[ends[1]] };
}

string poem::nextFour(string s, int rhythmLenth, int rhythmRegularity) {
	updateRythm(rhythmLenth, rhythmRegularity);
	updateText(s);
	if (nextRythmedPatternsPQ.size()>1) {
		prevRythmedFours.push(nextRythmedPatternsPQ.top());
		nextRythmedPatternsPQ.pop();
	}
	else {
		int prevSize = prevRythmedFours.size() + nextRythmedPatternsPQ.size();
		string res;
		do {
			res = nextWord();
		} while (nextRythmedPatternsPQ.size()==prevSize);
		for (int i = 0; i < prevSize; i++) {
			prevRythmedFours.push(nextRythmedPatternsPQ.top());
			nextRythmedPatternsPQ.pop();
		}
	}
	return bestSentenceOfBestFour();
}

string poem::prevFour() {
	if (prevRythmedFours.size()) {
		nextRythmedPatternsPQ.push(prevRythmedFours.top());
		prevRythmedFours.pop();
	}
	return bestSentenceOfBestFour();
}

string poem::nextPattern(int lineToChange) {
	currentPattern.nextPattern(lineToChange);
	return bestSentenseForPattern();
}

string poem::prevPattern(int lineToChange) {
	currentPattern.prevPattern(lineToChange);
	return bestSentenseForPattern();
}

void patternType::nextPattern(int lineIndex) {
	if (nextRythmedPatternsPQ[lineIndex].size()>1) {
		prevRythmedPatterns[lineIndex].push(nextRythmedPatternsPQ[lineIndex].top());
		nextRythmedPatternsPQ[lineIndex].pop();
	}
}

void patternType::prevPattern(int lineIndex) {
	if (prevRythmedPatterns[lineIndex].size()) {
		nextRythmedPatternsPQ[lineIndex].push(prevRythmedPatterns[lineIndex].top());
		prevRythmedPatterns[lineIndex].pop();
	}
}

string poem::bestSentence() {
	for (int slot = 0; slot < rythmPattern.size(); slot++) {
		int rp = rythmPattern[slot];
		addedWords[slot] = orderedWords[syllToWordsOrdered[slot][rp][0]] % outputSampleSize; // preliminary fill the array based on sylls
	}
	for (int i = 0; i < 4; i++) { // preliminary fill the array based on sylls/ends
		int slot = slots[i+1];
		int rp = rythmPattern[slot];
		int end = ends[i & 1];
		addedWords[slot] = orderedWords[syllEndToWordsOrdered[slot][rp][end][0]] % outputSampleSize; // preliminary fill the array based on sylls
	}
	for (int i = rythmPattern.size(); i < syllAndEndToWords.size(); i++)
		addedWords[i]=indexOfSpace;
	return result();
}

/*string poem::nextOption() {
	// look through all possible combinations of indices not exceeding the worse index:
	// e.g. slot1: 0,1,2, slot2: 3,4,5, slot3: 6,7,8. First combination is 0,0,0. The worst is 6 in slot3.
	// We need to look through all combinations of slot1 and slot2: 00,01,02,10,11,12,20,21,22.
	// No need to update pattern/slotOfAddedWord here.
	// Then - update the worse: 001. And again look through all combinations, etc.
	// When updating the worse - look through all the patterns/insert positions.

	// Looking through combinations - recousively (to move from the best combination to the worst):
	//    choose the best next word from all current slots
	//    if size of current slots is more than one - remove this slot from current slot and recourse
	//    else - move through all words from the last slot not exceeding weight from the previous slot

	// There could easily be a vast number of options. It is probably doen't make sense to look through every option.
	// It is more efficient to look through words in every slot (sorted).
	return"";
}*/

bool poem::rythmed(int rythmedBefore, int ap, int sn) {
	int rb = rythmedBefore % syllableInfo.total;
	if (rb + sn > syllableInfo.total)
		return false;
	if (sn < 2)
		return true;
	int accentPos = (rb + syllableInfo.regularity - syllableInfo.firstAccent+ap) % syllableInfo.regularity; // accentPos for words before
	return 0 == accentPos;
}

void poem::updatePermittedEnds(int slot0, int slot1, int syll1, map<int, set<int>> &permittedEnds, int slot, map<int, array<int, 3>> *scenTables) {
	permittedEnds[slot1].insert(syll1);
	// reset start slots for dot scenarios if slot<endSlot
	if (slot>=slot0 && slot <= slot1) {
		for (auto itScen = scenTables->begin(); itScen != scenTables->end(); itScen++) {
			auto scenTable = itScen->second;
			for (int i = 0; i < 3; i++) {
				if (scenTable[i] != -2) {
					scenTable[i] = -1;
				}
			}
		}
	}
}

string poem::excludeWord(int sentence, int slot) {
	//assert(sentence == 0);
	//assert(accumulatedLengths[17][6][1].size());
	if (slot == -1)
		return"";
	slot = slots[sentence] + 1 + slot;
	int word = addedWords[slot];
	if (slot & 1 && word == nodes[slot / 2].index){
		return bestSentenceOfBestFour();
	}
	// remove from syllAndEndToWords
	int syll = indexToSlots[word].rythmIndex;
	int end = indexToSlots[word].endIndex;
	auto itSyll = syllAndEndToWords[slot].find(syll);
	itSyll->second[end].erase(word);

	// remove index of word from syllToWordsOrdered, syllEndToWordsOrdered
	// it is difficult to remove from vector, but it is still the best choice
	// how to find index of word ? unordered_map is the best
	int index = wordToOrderedIndex[word + outputSampleSize * slot];
	auto it = find(syllToWordsOrdered[slot][syll].begin(), syllToWordsOrdered[slot][syll].end(), index);
	syllToWordsOrdered[slot][syll].erase(it);
	auto it2 = find(syllEndToWordsOrdered[slot][syll][end].begin(), syllEndToWordsOrdered[slot][syll][end].end(), index);
	syllEndToWordsOrdered[slot][syll][end].erase(it2);

	if (!rythms[syll][1]) {
		// update zeroCurr, zeroPrev, ZeroPrevPrev
		for (auto syll : syllAndEndToWords[slot]) {
			if (!rythms[syll.first][1]) {
				zeroCurr[slot] = false;
				zeroPrev[slot + 1] = false;
				zeroPrevPrev[slot + 2] = false;
				break;
			}
		}

		if (word == indexOfSpace) {
			// recalculate nextDotPosPossible, prevDotPosPossible, spacePermitted
			spacePermitted[slot] = false;
			int ndp = prevDotPos[slot];// update nextDotPos
			for (int i = slot; i > ndp; i--) {
				nextDotPosPossible[i] = nextDotPosPossible.size();
			}
			/*ndp = nextDotPos[slot];
			for (int i = slot; i < ndp; i++) {
				prevDotPosPossible[i] = -1;
			}*/
		}
	}

	if (!itSyll->second[end].size()) { // if no more words
		itSyll->second.erase(end); // erasing end could result in impossibility of some fours - will be checked when deleting pairs
		// check that end scenario 1 is possible
		bool endScenarioAffected = syllAndEndToWords[slot][syll].size() <= (slot & 1);
		if (!itSyll->second.size()) {// remove cascade
			syllAndEndToWords[slot].erase(itSyll);
			//syllToWordsOrdered[slot].erase(syll);
			//syllEndToWordsOrdered[slot].erase(syll);
		}
		// delete pairs linked to the end
		for (auto itSl0 = pairs.begin(); itSl0 != pairs.end(); ) {
			for (auto itSl1 = itSl0->second.begin(); itSl1 != itSl0->second.end();) {
				for (auto itSyll0 = itSl1->second.begin(); itSyll0 != itSl1->second.end(); ) {
					for (auto itSyll1 = itSyll0->second.begin(); itSyll1 != itSyll0->second.end(); ) {
						if (itSl0->first == slot && itSyll0->first == syll || itSl1->first == slot && itSyll1->first == syll) {
							itSyll1->second.erase(end); // affects standard/non-standard end scenario
						}
						if (!itSyll1->second.size()) { // affects fours
							auto it = next(itSyll1);
							itSyll0->second.erase(itSyll1);
							itSyll1 = it;
							continue;
						}
						itSyll1++;
					}
					if (!itSyll0->second.size()) {
						auto it = next(itSyll0);
						itSl1->second.erase(itSyll0);
						itSyll0 = it;
						continue;
					}
					itSyll0++;
				}
				if (!itSl1->second.size()) {
					auto it = next(itSl1);
					itSl0->second.erase(itSl1);
					itSl1 = it;
					continue;
				}
				itSl1++;
			}
			if (!itSl0->second.size()) {
				auto it = next(itSl0);
				pairs.erase(itSl0);
				itSl0 = it;
				continue;
			}
			itSl0++;
		}
		// delete from endToSlotToSyll
		endToSlotToSyll[end][slot].erase(syll);
		if (!endToSlotToSyll[end][slot].size())
			endToSlotToSyll[end].erase(slot);
		if (!endToSlotToSyll[end].size())
			endToSlotToSyll.erase(end);

		// Delete from fours 
		// Deleting end may result in some line patterns are not accessible. Is it possible to delete from fours ?
		// Start deleting from 0 slot. For 1 slot check that every start slot/syll exist in the end of 0 slot.
		// Check whether the end scenario is possible. Check that the both pairs exist
		// Would it be better to create fours from zero ? It is worse - every pair should be compared against every other
		map<int, set<int>> permittedStarts;// = { 0, {0} };
		permittedStarts[-1].insert(0);
		map<int, set<int>> permittedEnds;
		map<int, set<int>> permittedPairEnds[2];
		map<int, set<int>>::iterator itPairEnds;// for 2 and 3 - check that slot/syll is within a set of pair ends corresponding to 0 and 1 slot/sylls
		map<int, set<int>> permittedPairStarts;
		map<int, set<int>>::iterator itPairs;// for 0 and 1 - check that slot/syll is within set of pair starts
		// create permittedPairStarts
		for (auto sl0 : pairs) {
			for (auto sl1 : sl0.second) {
				for (auto syll0 : sl1.second) {
					permittedPairStarts[sl0.first].insert(syll0.first);
				}
			}
		}
		for (int i = 0; i < 4; i++) {
			for (auto itSl0 = fours[i].begin(); itSl0 != fours[i].end(); ) {
				// check that slot0/syll0 is included into ends of the previous line
				auto itSl0_permitted = permittedStarts.find(itSl0->first-1);
				if (itSl0_permitted != permittedStarts.end()) {
					for (auto itSl1 = itSl0->second.begin(); itSl1 != itSl0->second.end(); ) {
						// check that the pair exists: if i<2 - check start of pair, else - pair in full
						if (i < 2) { // check that slot1/syll1 is included in pairs start
							itPairs = permittedPairStarts.find(itSl1->first);
							if (itPairs == permittedPairStarts.end()) {
								itSl1++;
								continue;
							}
						}
						else { // check that slot1/syll1 is included in permittedPairEnds
							itPairEnds = permittedPairEnds[i & 1].find(itSl1->first);
							if (itPairEnds == permittedPairEnds[i & 1].end()) {
								itSl1++;
								continue;
							}
						}
						for (auto itSyll1 = itSl1->second.begin(); itSyll1 != itSl1->second.end();) {
							// check that the pair exists: if i<2 - check start of pair, else - pair in full
							if (i < 2) {
								auto itSyll = itPairs->second.find(itSyll1->first);
								if (itSyll == itPairs->second.end()) {// if there is no pair
									auto it = next(itSyll1);
									itSl1->second.erase(itSyll1); // erase from four
									itSyll1 = it;
									continue;
								}
							}
							else {
								auto itSyll = itPairEnds->second.find(itSyll1->first);
								if (itSyll == itPairEnds->second.end()) {// if there is no pair 
									auto it = next(itSyll1);
									itSl1->second.erase(itSyll1); // erase from four
									itSyll1 = it;
									continue;
								}
							}
							if (itSl1->first == slot && itSyll1->first == syll && endScenarioAffected) {
								itSyll1->second.erase(slot & 1);//erase
								if (itSyll1->second.size()) {// check that the pairs exist and add to permitted ends
									updatePermittedEnds(itSl0->first, itSl1->first, itSyll1->first, permittedEnds, slot, &itSyll1->second);
								}
								else {
									auto it = next(itSyll1);
									itSl1->second.erase(itSyll1); // erase from four
									itSyll1 = it;
									continue;
								}
							}
							else {// check that the pairs exist and add to permitted ends
								updatePermittedEnds(itSl0->first, itSl1->first, itSyll1->first, permittedEnds, slot, &itSyll1->second);
							}
							itSyll1++;
						}
						if (!itSl1->second.size()) {
							auto it = next(itSl1);
							itSl0->second.erase(itSl1);
							itSl1 = it;
							continue;
						}
						itSl1++;
					}
				}
				else { // erase slot
					auto it = next(itSl0);
					fours[i].erase(itSl0);
					itSl0 = it;
					continue;
				}
				if (!itSl0->second.size()) {
					auto it = next(itSl0);
					fours[i].erase(itSl0);
					itSl0 = it;
					continue;
				}
				itSl0++;
			}
			permittedStarts = permittedEnds;
			permittedEnds.clear();
			if (i < 2) {//fill in permittedPairEnds
				for (auto startSlot : permittedStarts) {
					for (auto startSyll : startSlot.second) {
						for (auto endSlot : pairs[startSlot.first]) {// incorrect: pair may not exist now
							for (auto endSyll : endSlot.second[startSyll]) {
								permittedPairEnds[i][endSlot.first].insert(endSyll.first);
							}
						}
					}
				}
			}
		}
		/*if (itSyll->second.size() == 0) {
			auto it = next(itSyll);
			syllAndEndToWords[slot].erase(syll); 
			itSyll = it;
			// if no syll - it results in impossibility of some patterns and could result in impossibility of some fours. All accLen tables should be recalculated

			if (itSyll == syllAndEndToWords[slot].end()) {
				// results in impossibility of all patterns
			}
		}*/
		// recreate all accLen tables - just shrink all existing tables up to the slot and updateAccumLengthsForNewEnd
		for (auto it = accumulatedLengths.begin(); it != accumulatedLengths.end(); it++) {
			int startSlot = it->first;
			int endSlot = it->second.size()+startSlot-1;
			if (slot>=startSlot && endSlot >= slot) {
				it->second.resize(slot - startSlot);
				updateAccumLengthsForNewEnd(&it->second, startSlot, endSlot);
			}
		}

		nextRythmedPatternsPQ = priority_queue<patternType>();
		prevRythmedFours = stack<patternType>();
		updateFourForNewSyllRecoursively(0, 0, true);
	}

	//assert(accumulatedLengths[17][6][1].size());
	if (nextRythmedPatternsPQ.size()) {
		return bestSentenceOfBestFour();
	}
	else {
		return nextWord();
	}
}

#include <iostream>

void poem::print2() {
	for (int i = 0; i < addedWords.size(); i++) {
		cout << indexToSlots[addedWords[i]].s << " ";
		cout << nodes[i].s << " ";
	}
	cout << "\n";
}

string poem::result() {
	string s;
	int lineIndex = 1;
	for (int i = 0; i < rythmPattern.size(); i++) {
		if (addedWords[i] != indexOf$)
			s += indexToSlots[addedWords[i]].s + " ";
		else
			s += nodes[i / 2].s + " ";
		if (i == slots[lineIndex]) {
			lineIndex++;
			if (i == rythmPattern.size() - 1) {
				for (int i = rythmPattern.size(); i < addedWords.size(); i++) {
					s += "~ ";
				}
			}
			s += "\r\n";
		}
	}
	return s;
}

vector<string> poem::possibleWords(int slot, int line) {
	vector<string> result;
	int s = slots[line]+slot+1;
	if (s < syllAndEndToWords.size()) {
		int rp;
		if (s < rythmPattern.size())
			rp = rythmPattern[s];
		else
			rp = rythms[0][0];
		if (slots[line + 1] == s) {
			int end = ends[line & 1];
			for (int index = 0; index < syllEndToWordsOrdered[s][rp][end].size(); index++) {
				result.push_back(indexToSlots[orderedWords[syllEndToWordsOrdered[s][rp][end][index]] % outputSampleSize].s);
			}
		}
		else {
			for (int index = 0; index < syllToWordsOrdered[s][rp].size(); index++) {
				result.push_back(indexToSlots[orderedWords[syllToWordsOrdered[s][rp][index]] % outputSampleSize].s);
			}
		}
	}
	return result;
}

string poem::updateEnd(int endIndex, string end) {
	ends[endIndex] = nodeType::endStringToCode[end];
	return bestSentence();
}