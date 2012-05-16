/*
 Copyright (c) 2004 Markus Niemist�

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or
 sell copies of the Software, and to permit persons to whom
 the Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall
 be included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
*/

#ifndef __TABLE_HH
#define __TABLE_HH

#include <Fl/Fl.H>
#include <Fl/Fl_Box.H>
#include <Fl/Fl_Group.H>
#include <Fl/Fl_Menu_Item.H>
#include <Fl/Fl_Scrollbar.H>
#include <Fl/Fl_Widget.H>
#include <vector>

#define TABLE_WHEN_DCLICK		16


class Table : public Fl_Group {
public:
	typedef bool (*Highlighter)(int, const char **, Fl_Color &);
	typedef int (*Comparator)(const char*, const char*);

private:
	struct ColumnInfo {
		bool hidden;
		const char *title;
		int width;
		Fl_Align align;
		Comparator comparator;
	};

	// Scrollbars
	Fl_Scrollbar *hScroll, *vScroll;

	// Popup menu
	const Fl_Menu_Item *popupMenu;
	bool menuAlloc;

	// Column data
	std::vector<struct ColumnInfo> header;

	// Cell data
	std::vector<char**> data;
	Highlighter highlighter;

	// Table dimensions
	int tableHeight, tableWidth;
	int oX, oY, oW, oH;	// Outer dimensions (widget - box)
	int iX, iY, iW, iH;	/*
				 * Table area dimensions
				 * (outer dimension - header - scrollbars)
				 */

	// For optimization
	int topRow, bottomRow, leftCol, rightCol;
	int topRowY, leftColX;

	int nCols, nRows;	// Number of rows and columns
	int cPos;		// Column where new entry is added.

	int resizing, dragX;
	int pushed;
	int sortColumn;

	// Object sizes
	int scrollbarSize;
	int headerHeight;
	int rowHeight;

	int selected;
	char **curRow;

	// Various flags
	bool ascent;
	bool canResize, canSort;
	bool noMoreColumns;
	bool toBeSorted;
	bool dimensionsChanged;
	bool headerEnabled;

	void dSort(int start, int end, Comparator compare);
	void aSort(int start, int end, Comparator compare);

protected:
	virtual int handle(int event);

	virtual void drawHeader(int x, int y);
	virtual void drawRow(int row, char *rowData[], int x, int y);

	virtual void draw();
	virtual void resize(int x, int y, int w, int h);

	void calcDimensions();
	void scrolled();
	void resized();

	static void scrollCallback(Fl_Widget *widget, void *data);

public:
	Table(int x, int y, int w, int h, char *label = NULL);
	~Table();

	bool headerOn() const;
	void headerOn(bool enabled);
	bool allowResize() const;
	void allowResize(bool allow);
	bool allowSort() const;
	void allowSort(bool allow);

	int headerSize() const;
	void headerSize(int height);
	int rowSize() const;
	void rowSize(int height);
	int scrollbSize() const;
	void scrollbSize(int size);

	Fl_Align columnAlign(int column) const;
	void columnAlign(int column, Fl_Align align);
	int columnWidth(int column) const;
	void columnWidth(int column, int width);
	const char *columnTitle(int column);
	void columnTitle(int column, const char *title);
	bool columnHidden(int column);
	void columnHidden(int column, bool hidden);

	void sort();
	void sort(int column, bool ascent);
	void getSort(int &sortColumn, bool &ascent);

	void setHighlighter(Highlighter highlighter);

	void addColumn(const char *label, int width = 150,
	    Fl_Align align = (Fl_Align)(FL_ALIGN_LEFT | FL_ALIGN_CLIP),
	    Comparator comparator = NULL);
	void addHiddenColumn(const char *label);

	void addCell(const char *data);
	void addCell(int data);
	void addRow(int cols, ...);
	void addFromTSV(const char *data);
	void removeRow(int row);
	void clear(bool removeColumns = false);

	void where(int x, int y, int &row, int &column, int &resize);
	void scrollTo(int pos);

	int columns();
	int rows();
	void value(int selection);
	int value();
	char *valueAt(int row, int column);
	int intValueAt(int row, int column);
	void valueAt(int row, int column, char *data);
	void valueAt(int row, int column, int data);

	const char **getRow(int row);

	const Fl_Menu_Item *menu();
	void menu(const Fl_Menu_Item *m);
	void menuCopy(const Fl_Menu_Item *m);
	void menuClear();

	static int compareInt(const char *val1, const char *val2);
};

#endif
