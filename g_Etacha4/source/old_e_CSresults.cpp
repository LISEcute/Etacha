//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#include <vcl.h>
#pragma hdrstop

#include "e_CSresults.h"
#include "Constants.h"
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#pragma package(smart_init)
#pragma resource "*.dfm"
#include "e_myextern.h"
#include <math.h>
#include <string.h>

TFormResults *FormResults=0;

extern char* ElementName(int IZ);
extern double gAb,gZb,gEnergy;
extern double gZt;


const int I2_Maska[6][6] = {
//2s,2p,3s,3p,3d,4n
 {15,16,17,18,19,20}, //1s
 { 0,32,21,22,23,24}, //2s
 { 0, 0,25,26,27,28}, //2p
 { 0, 0, 0,33, 0,29}, //3s
 { 0, 0, 0, 0,34,30}, //3p
 { 0, 0, 0, 0, 0,31}};//3d

 const int E2_Maska[7][7] = {
//2s,2p,3s,3p,3d,4n,5n
 {29,30,31,32,33,34,22}, //1s
 { 0,46,35,36,37,38,23}, //2s
 { 0, 0,39,40,41,42,24}, //2p
 { 0, 0, 0,47, 0,43,25}, //3s
 { 0, 0, 0, 0,48,44,26}, //3p
 { 0, 0, 0, 0, 0,45,27}, //3d
 { 0, 0, 0, 0, 0, 0,28}  //n4
 };

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
__fastcall TFormResults::TFormResults(TComponent* Owner)
        : TForm(Owner)
{
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::FormCreate(TObject *Sender)
{
        if(EtachaVersion >= etacha_v45)
                        {
                        StringGrid_E2->RowCount = 8;              // right bottom
                        StringGrid_E2->ColCount = 8;
                        StringGrid_E2->DefaultRowHeight = 20;

                        StringGrid_I2->RowCount = 7;              // right top
                        StringGrid_I2->ColCount = 7;
                        StringGrid_I2->DefaultRowHeight = 23;
                        }
        else if(EtachaVersion <= etacha_v3)
                        {
                        StringGrid_E2->RowCount = 6;
                        StringGrid_E2->ColCount = 6;
                        StringGrid_E2->DefaultRowHeight = 27;

                        StringGrid_I2->RowCount = 6;
                        StringGrid_I2->ColCount = 6;
                        StringGrid_I2->DefaultRowHeight = 27;
                        }
        else            {
                        StringGrid_E2->RowCount = 7;
                        StringGrid_E2->ColCount = 7;
                        StringGrid_E2->DefaultRowHeight = 23;

                        StringGrid_I2->RowCount = 7;
                        StringGrid_I2->ColCount = 7;
                        StringGrid_I2->DefaultRowHeight = 23;
                        }

      //-----------------------------------------------------------------
        if(EtachaVersion >= etacha_v34)
                        {
                        StringGrid_Inform->RowCount = 8;                  // "inform"  top-left
                        StringGrid_Inform->DefaultRowHeight = 20;

                        StringGrid_EditCS->RowCount = 8;                 //  "EdiCS"  bottom left
                        StringGrid_EditCS->DefaultRowHeight = 20;
                        }
        else            {
                        StringGrid_Inform->RowCount = 7;
                        StringGrid_Inform->DefaultRowHeight = 23;

                        StringGrid_EditCS->RowCount = 7;
                        StringGrid_EditCS->DefaultRowHeight = 23;
                        }



StringGrid_EditCS->Color = TColor(RGB(255, 230, 204));
StringGrid_EditCS->FixedColor = TColor(RGB(255, 128, 0));

StringGrid_Inform->Color = TColor(RGB(250, 220, 214));
StringGrid_Inform->FixedColor = TColor(RGB(250, 118, 20));

StringGrid_I2->Color = TColor(RGB(245, 210, 224));
StringGrid_I2->FixedColor = TColor(RGB(245, 108, 40));

StringGrid_E2->Color = TColor(RGB(240, 200, 234));
StringGrid_E2->FixedColor = TColor(RGB(240, 98, 60));

const int   ColW_CS[] = {30,53,90,90,90};
const char *ColN_CS[] = {"  N","(sub)shell","MEC (capture)","REC (capture)","IONization"};
const char *Shells7[] = {"1s","2s","2p","3s","3p","3d","n=4","n=5,6"};

const int   Col2W_CS[] = {30,52,130,140};
const char *Col2N_CS[] = {"  N","(sub)shell","Capture (MEC+REC)","IONization (includes n>4)"};

const int   Col2IW_CS[] = {80,80,80,80,80,80,80,80};
const char *Versions[] = {"23","3","34","4","45"};

const char *chIon[] = {"CDW-EIS","PWBA"};
const char *chExc[] = {"SE","PWBA"};
int total = 0;
AnsiString Buffer = Versions[EtachaVersion];
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWwww   Labels
LabelVersion->Caption = "Version " + Buffer;
Label_Ionization->Caption =chIon[IonizationModel];
Label_Excitation->Caption =chExc[ExcitationModel];

char Cbeam[3];
char Ctarg[3];
strncpy(Cbeam,ElementName(gZb),2); int kb=2; if(Cbeam[1]==' ')kb=1; Cbeam[kb]=0;
strncpy(Ctarg,ElementName(gZt),2); int kt=2; if(Ctarg[1]==' ')kt=1; Ctarg[kt]=0;

Buffer.sprintf("%.0f%s (%.1fMeV/u) + %s",gAb,Cbeam, gEnergy,Ctarg);
Label_Reaction->Caption = Buffer;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW   Edit 1
for(int i=0;i<StringGrid_EditCS->ColCount; i++)
        {
        StringGrid_EditCS->ColWidths[i]=ColW_CS[i];
        total += ColW_CS[i]+2;
        Buffer = ColN_CS[i];
        if(i>=2) Buffer = "    " + Buffer;
        StringGrid_EditCS->Cells[i][0] =  Buffer;
        }
StringGrid_EditCS->Width=total;
StringGrid_EditCS->CellRect(2,2);

for(int i=0;i<StringGrid_EditCS->RowCount; i++)
        {
        int row = i+1;
        StringGrid_EditCS->Cells[0][row] = row;
        Buffer = Shells7[i];   Buffer = "    " + Buffer;
        StringGrid_EditCS->Cells[1][row]=  Buffer;
        }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW   inform 1
total = 0;
for(int i=0;i<StringGrid_Inform->ColCount; i++)
        {
        StringGrid_Inform->ColWidths[i]=Col2W_CS[i];
        total += Col2W_CS[i]+2;
        Buffer = Col2N_CS[i];
        if(i>=2) Buffer = "    " + Buffer;
        StringGrid_Inform->Cells[i][0] =  Buffer;
        }
StringGrid_Inform->Width=total;

for(int i=0;i<StringGrid_Inform->RowCount; i++)
        {
        int row = i+1;
        StringGrid_Inform->Cells[0][row] = row;
        Buffer = Shells7[i];   Buffer = "   " + Buffer;
        StringGrid_Inform->Cells[1][row]=  Buffer;
        }
StringGrid_Inform->CellRect(2,2);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW   inform 2
total = 0;
#define ColWidth_l2 70
for(int i=0;i<StringGrid_I2->ColCount; i++)
        {
        StringGrid_I2->ColWidths[i]=ColWidth_l2;
        total += ColWidth_l2+2;

        if(i==0)  StringGrid_I2->Cells[i][0] =  "  From / To";
        else    {
                Buffer = Shells7[i];
                Buffer = "       " + Buffer;
                StringGrid_I2->Cells[i][0] =  Buffer;
                }
        }
StringGrid_I2->Width=total-2;

for(int i=0;i<StringGrid_I2->RowCount; i++)
        {
        int row = i+1;
        Buffer = Shells7[i];   Buffer = "       " + Buffer;
        StringGrid_I2->Cells[0][row]=  Buffer;
        }
StringGrid_I2->CellRect(1,1);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW   Edit2 2
total = 0;
#define ColWidth_E2 70
for(int i=0;i<StringGrid_E2->ColCount; i++)
        {
        StringGrid_E2->ColWidths[i]=ColWidth_E2;
        total += ColWidth_E2+2;

        if(i==0)  StringGrid_E2->Cells[i][0] =  "  From / To";
        else    {
                Buffer = Shells7[i];
                Buffer = "       " + Buffer;
                StringGrid_E2->Cells[i][0] =  Buffer;
                }
        }
StringGrid_E2->Width=total-3;

for(int i=0;i<StringGrid_E2->RowCount; i++)
        {
        int row = i+1;
        Buffer = Shells7[i];   Buffer = "       " + Buffer;
        StringGrid_E2->Cells[0][row]=  Buffer;
        }
StringGrid_E2->CellRect(1,1);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW Values
//--------------------------------------------------------   inform 1

for(int row=0; row<7; row++)
   for(int col=0; col<2; col++)
               StringGrid_Inform->Cells[col+2][row+1] =
                   Buffer.FormatFloat("     0.0000e+0",GSEC[col*7+row+1]);


//--------------------------------------------------------   inform 2

for(int row=0; row<6; row++)
   for(int col=0; col<6; col++)
        {
        int Flag = I2_Maska[row][col];
        if(Flag > 0)
               StringGrid_I2->Cells[col+1][row+1] = Buffer.FormatFloat("  0.0000e+0",GSEC[Flag]);
        }

//--------------------------------------------------------   edit 1

for(int row=0; row<7; row++)
   for(int col=0; col<3; col++)
               StringGrid_EditCS->Cells[col+2][row+1] =
                   Buffer.FormatFloat("     0.0000e+0",Gsecs[col*7+row+1]);

//--------------------------------------------------------   edi2 2

for(int row=0; row<7; row++)
   for(int col=0; col<7; col++)
        {
        int Flag = E2_Maska[row][col];
        if(Flag > 0)
               StringGrid_E2->Cells[col+1][row+1] = Buffer.FormatFloat("  0.0000e+0",Gsecs[Flag]);
        }

BtnResults_Accept->Enabled = false;

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::StringGrid_I2DrawCell(TObject *Sender,
      int ACol, int ARow, TRect &Rect, TGridDrawState State)
{
if(ACol <= 0 || ARow <= 0) return;
//[column][row]
const int ForbidMaska[6][6] = {
{1,1,1,1,1,1},
{0,1,1,1,1,1},
{0,0,1,1,1,1},
{0,0,0,1,0,1},
{0,0,0,0,1,1},
{0,0,0,0,0,1}};

int lRow = ARow-1;
int lCol = ACol-1;
int Flag = ForbidMaska[lRow][lCol];

if(!Flag)
        {
        StringGrid_I2->Canvas->Brush->Color = clLtGray;
        StringGrid_I2->Canvas->FillRect(Rect);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::StringGrid_E2DrawCell(TObject *Sender,
      int ACol, int ARow, TRect &Rect, TGridDrawState State)
{
if(ACol <= 0 || ARow <= 0) return;

const int ForbidMaskaE2[7][7] = {
{1,1,1,1,1,1,1},
{0,1,1,1,1,1,1},
{0,0,1,1,1,1,1},
{0,0,0,1,0,1,1},
{0,0,0,0,1,1,1},
{0,0,0,0,0,1,1},
{0,0,0,0,0,0,1}};


int lRow = ARow-1;
int lCol = ACol-1;
int Flag = ForbidMaskaE2[lRow][lCol];

if(!Flag)
        {
        StringGrid_E2->Canvas->Brush->Color = clLtGray;
        StringGrid_E2->Canvas->FillRect(Rect);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::BtnResults_AcceptClick(TObject *Sender)
{
        //--------------------------------------------------------   edit 1

for(int row=0; row<7; row++)
   for(int col=0; col<3; col++)
               {
               int index = col*7+row+1;
               double v = fabs(atof(StringGrid_EditCS->Cells[col+2][row+1].c_str()));
               Gcor[index]= Gseci[index] > 0 ? v/Gseci[index] : 1;
               }


        //--------------------------------------------------------   edi2 2

for(int row=0; row<7; row++)
   for(int col=0; col<7; col++)
        {
        int Flag = E2_Maska[row][col];
        if(Flag > 0)
              {
              double v = fabs(atof(StringGrid_E2->Cells[col+1][row+1].c_str()));
              Gcor[Flag]= Gseci[Flag] > 0 ? v/Gseci[Flag] : 1;
              }
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::StringGrid_E2GetEditText(TObject *,
      int , int , AnsiString &)
{
BtnResults_Accept->Enabled = true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void __fastcall TFormResults::StringGrid_EditCSGetEditText(TObject *,
      int , int , AnsiString &)
{
BtnResults_Accept->Enabled = true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

