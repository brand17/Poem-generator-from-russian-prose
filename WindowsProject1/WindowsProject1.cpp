#define NOMINMAX
#include <windows.h>
#include "poem.h"
#include <string>
#include "resource.h"
poem* p;

#include <commctrl.h>
using namespace std;

typedef struct SubButtonData_tag SubButtonData;
struct SubButtonData_tag {
	// ... data for the subclass implementation
};

struct coords {
	int x;
	int y;
	string word;
};

coords GetEditWord(HWND hwndEdit, POINT point)
{
	int nChar = (int)::SendMessage(hwndEdit, EM_CHARFROMPOS, 0, MAKELPARAM(point.x, point.y));
	int nLineIndex = HIWORD(nChar);
	int nCharIndex = LOWORD(nChar);

	int nLineStart = (int)::SendMessage(hwndEdit, EM_LINEINDEX, nLineIndex, 0);
	// тут длина строки текста 
	int nLineLength = (int)::SendMessage(hwndEdit, EM_LINELENGTH, nLineStart, 0);
	int nLineEnd = nLineStart + nLineLength; // конец строки
	if (nCharIndex >= nLineStart + nLineLength || nCharIndex < 0) return { -1,nLineIndex,"" };

	// получаем текст строки
	char* cc = new char[nLineLength];
	*(LPWORD)cc = (WORD)nLineLength;
	::SendMessage(hwndEdit, EM_GETLINE, nLineIndex, (LPARAM)cc);

	int nBufChar = nLineLength - (nLineEnd - nCharIndex);
	int res = 0;
	for (int i = nBufChar; i; i--) {
		if (cc[i] == char(' ')) res++;
	}
	int first = nBufChar;
	do {
		if (cc[first] != char(' ')) first--;
		else {
			first++;
			break;
		}
	} while (first > 0);
	int end = first;
	do {
		if (cc[end] != char(' ')) end++;
		else break;
	} while (end < nLineLength);
	string s = string(&cc[first], end - first);
	delete[] cc;

	return { res, nLineIndex, s };// sReturn;
}

// The API treats the pair (uSubclassId, SubButtonProc) as unique identification
// of the subclass. Assuming we do not need multiple subclass levels of
// the same control (which would share the subclass procedure), we do not need
// to deal much with the ID. Only if we would need such esoteric generality, we
// would use it for distinguishing among the subclasses.
static UINT_PTR uFinalEditViewId = 0;
static UINT_PTR uBaseEditViewId = 1;
coords wordToDelete;
#define ID_FILE_EXIT 9001

void updateEndComboBoxes(HWND hwnd) {
	auto ends = p->currentEnds();
	int index = SendDlgItemMessage(hwnd, IDC_COMBO1, CB_RESETCONTENT, 0, 0);
	index = SendDlgItemMessage(hwnd, IDC_COMBO2, CB_RESETCONTENT, 0, 0);
	int combos[] = { IDC_COMBO1 ,IDC_COMBO2 };
	int endIndex[2] = { 0 };
	for (int i = 0; i < 2; i++) {
		for (auto e : p->possibleEnds(i)) {
			int index = SendDlgItemMessage(hwnd, combos[i], CB_ADDSTRING, 0, (LPARAM)e.c_str());
			if (ends[i] == e)
				SendDlgItemMessage(hwnd, combos[i], CB_SETCURSEL, (WPARAM)endIndex[i], (LPARAM)0);
			endIndex[i]++;
		}
	}
}

static LRESULT CALLBACK
finalEditViewProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam, UINT_PTR uSubclassId, DWORD_PTR dwData)
{
	SubButtonData* lpData = (SubButtonData*)dwData;

	switch (uMsg) {
		// Handle all messages we want to customize:
		//...
	case WM_RBUTTONDOWN: {
		POINT scr;
		GetCursorPos(&scr);
		POINT loc = scr;
		ScreenToClient(hwnd, &loc);
		wordToDelete = GetEditWord(hwnd, loc);
		HMENU hSubMenu = CreatePopupMenu();
		LPSTR cString = _strdup(wordToDelete.word.c_str());
		AppendMenu(hSubMenu, MF_STRING, ID_DELETE_WORD_MENU_ITEM, cString);
		TrackPopupMenu(hSubMenu, TPM_BOTTOMALIGN | TPM_LEFTALIGN, scr.x, scr.y, 0, hwnd, NULL);
		DestroyMenu(hSubMenu);
		return 0;
	}
	case WM_LBUTTONDOWN: {
		POINT scr;
		GetCursorPos(&scr);
		POINT loc = scr;
		ScreenToClient(hwnd, &loc);
		wordToDelete = GetEditWord(hwnd, loc);
		if (wordToDelete.x == -1)
			return 0;
		vector<string> words = p->possibleWords(wordToDelete.x, wordToDelete.y);
		HMENU hSubMenu = CreatePopupMenu();
		for (int i = 0; i < words.size(); i++) {
			LPSTR cString = _strdup(words[i].c_str());
			AppendMenu(hSubMenu, MF_STRING, ID_FILE_EXIT, cString);
		}
		TrackPopupMenu(hSubMenu, TPM_BOTTOMALIGN | TPM_LEFTALIGN, scr.x, scr.y, 0, hwnd, NULL);
		DestroyMenu(hSubMenu);

		return 0;
	}
	case WM_COMMAND:
		switch (LOWORD(wParam)) {
		case ID_DELETE_WORD_MENU_ITEM:
			string result = p->excludeWord(wordToDelete.y, wordToDelete.x);
			updateEndComboBoxes(GetParent(hwnd));
			SetDlgItemText(GetParent(hwnd), IDC_EDIT2, result.c_str());
		}
		break;
	}

	return DefSubclassProc(hwnd, uMsg, wParam, lParam);
}

INT_PTR CALLBACK DlgProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	int len;
	string result;
	switch (Message)
	{
	case WM_INITDIALOG:
		SetWindowSubclass(GetDlgItem(hwnd, IDC_EDIT2), finalEditViewProc, uFinalEditViewId, NULL);
		//SetDlgItemText(hwnd, IDC_EDIT1, "п`альцы сж`али . свет появ`ился . мысль пришл`а . судьб`а ждал`а .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "я всегд`а исп`ытываю сч`астье там . где я мог`у быть насто`ящим .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "но что бы ты ни предрек`ал, иисус назар`янин, напрасны усилия твои.");
		SetDlgItemText(hwnd, IDC_EDIT1, "к`аждый жив`ёт . как х`очет . и распл`ачивается за это сам .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "иногд`а мом`ент . кот`орый ты так д`олго ждал . прих`одит в с`амое неподход`ящее вр`емя .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "не соверш`ай класс`ическую ош`ибку всех `умников . не д`умай . что нет люд`ей умн`ее теб`я .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "новог`однее настро`ение . `это когд`а рад в`идеть д`аже тех . кто ош`ибся дв`ерью .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "челов`ека выда`ёт . то над чем сме`ётся он .");//short phrase - low probability of good poem
		//SetDlgItemText(hwnd, IDC_EDIT1, "будь счастл`ивым . в `этот с`амый миг . он и есть . жизнь тво`я .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "неуд`ача . `это пр`осто возм`ожность нач`ать сн`ова . но уж`е б`олее м`удро .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "пом`имо ч`удом уцел`евшего д`ома сыс`оева . в печ`атниковом и сос`еднем колок`ольниковом пере`улках б`ыло ещ`ё н`есколько не столь ц`енных . но . тем не м`енее . очаров`ательных дерев`янных особнячк`ов .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "Мои друзь`я . спас`ибо за т`ёплые слов`а . `очень мне при`ятно .");
		//SetDlgItemText(hwnd, IDC_EDIT1, "Я думаю о тебе так много, что мне даже странно, откуда берётся время на всё остальное.");
		//SetWindowSubclass(GetDlgItem(hwnd, IDC_EDIT1), baseEditViewProc, uBaseEditViewId, NULL);
		SetDlgItemText(hwnd, IDC_EDIT5, to_string(p->syllableInfo.total).c_str());
		SetDlgItemText(hwnd, IDC_EDIT6, to_string(p->syllableInfo.regularity).c_str());
		break;
	case WM_COMMAND:
		switch (LOWORD(wParam)) {
		case IDC_COMBO1: case IDC_COMBO2:
			if (HIWORD(wParam) == CBN_SELCHANGE)
				// If the user makes a selection from the list:
				//   Send CB_GETCURSEL message to get the index of the selected list item.
				//   Send CB_GETLBTEXT message to get the item.
				//   Display the item in a messagebox.
			{
				int ItemIndex = SendMessage((HWND)lParam, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
				TCHAR  ListItem[256];
				(TCHAR)SendMessage((HWND)lParam, (UINT)CB_GETLBTEXT,
					(WPARAM)ItemIndex, (LPARAM)ListItem);
				string s = p->updateEnd(LOWORD(wParam) == IDC_COMBO2, string(ListItem));
				SetDlgItemText(hwnd, IDC_EDIT2, s.c_str());
			}
			break;
		case IDC_BUTTON1: case IDC_BUTTON2:
			len = GetWindowTextLength(GetDlgItem(hwnd, IDC_EDIT1));
			if (len > 0)
			{
				char* buf;
				buf = (char*)GlobalAlloc(GPTR, len + 1);
				GetWindowTextA(GetDlgItem(hwnd, IDC_EDIT1), buf, len + 1);
				string result;
				HWND button = GetDlgItem(hwnd, LOWORD(wParam));
				EnableWindow(button, false);
				len = GetWindowTextLength(GetDlgItem(hwnd, IDC_EDIT5));
				char* bufLen = (char*)GlobalAlloc(GPTR, len + 1);
				GetWindowTextA(GetDlgItem(hwnd, IDC_EDIT5), bufLen, len + 1);
				len = GetWindowTextLength(GetDlgItem(hwnd, IDC_EDIT6));
				char* bufReg = (char*)GlobalAlloc(GPTR, len + 1);
				GetWindowTextA(GetDlgItem(hwnd, IDC_EDIT6), bufReg, len + 1);
				switch (LOWORD(wParam)) {
				case IDC_BUTTON1: // next pattern
					result = p->nextFour(buf, atoi(bufLen), atoi(bufReg));
					updateEndComboBoxes(hwnd);
					break;
				case IDC_BUTTON2: // previous pattern
					result = p->prevFour();
					if (result != "")
						updateEndComboBoxes(hwnd);
					break;
				}
				EnableWindow(button, true);
				if (result == "") {
					result = "No rythmes found";
				}
				SetDlgItemText(hwnd, IDC_EDIT2, result.c_str());
				GlobalFree((HANDLE)buf);
			}
			break;
		case IDC_BUTTON5:
			SetDlgItemText(hwnd, IDC_EDIT2, p->nextPattern(0).c_str());
			break;
		case IDC_BUTTON6:
			SetDlgItemText(hwnd, IDC_EDIT2, p->nextPattern(1).c_str());
			break;
		case IDC_BUTTON7:
			SetDlgItemText(hwnd, IDC_EDIT2, p->nextPattern(2).c_str());
			break;
		case IDC_BUTTON8:
			SetDlgItemText(hwnd, IDC_EDIT2, p->nextPattern(3).c_str());
			break;
		case IDC_BUTTON9:
			SetDlgItemText(hwnd, IDC_EDIT2, p->prevPattern(0).c_str());
			break;
		case IDC_BUTTON10:
			SetDlgItemText(hwnd, IDC_EDIT2, p->prevPattern(1).c_str());
			break;
		case IDC_BUTTON11:
			SetDlgItemText(hwnd, IDC_EDIT2, p->prevPattern(2).c_str());
			break;
		case IDC_BUTTON12:
			SetDlgItemText(hwnd, IDC_EDIT2, p->prevPattern(3).c_str());
			break;
		}
		break;
	case WM_CLOSE:
		EndDialog(hwnd, 0);
		break;
	default:
		return FALSE;
	}
	return TRUE;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	setlocale(LC_ALL, "");
	p = new poem();
	int r = DialogBox(hInstance, MAKEINTRESOURCE(IDD_DIALOG1), NULL, DlgProc);
	DWORD err = 0;
	if (r == -1) {
		err = GetLastError();
	}
	delete p;
	return true;
}
