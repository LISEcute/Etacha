#include <QTextEdit>
#include <QTextDocument>
#include <QTextCursor>
#include <QScrollBar>


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void scrollLogToBottom(QTextEdit *editWidget)
{
    QScrollBar* bar =  editWidget->verticalScrollBar();
    bar->setValue(bar->maximum());
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void fastAppend(QTextEdit *editWidget, const QString& message)
{
    const bool atBottom =  editWidget->verticalScrollBar()->value() ==
                           editWidget->verticalScrollBar()->maximum();
    QTextDocument* doc = editWidget->document();
    QTextCursor cursor(doc);
    cursor.movePosition(QTextCursor::End);
    cursor.beginEditBlock();
    cursor.insertBlock();
    cursor.insertHtml(message);
    cursor.endEditBlock();

    //scroll scrollarea to bottom if it was at bottom when we started
    //(we don't want to force scrolling to bottom if user is looking at a
    //higher position)
    if (atBottom) {
        scrollLogToBottom(editWidget);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
