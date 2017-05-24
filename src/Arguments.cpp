/****************************************************************************
**
** Copyright (C) 2011 Andrey Kartashov .
** All rights reserved.
** Contact: Andrey Kartashov (porter@porter.st)
**
** This file is part of the global module of the genome-tools.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General Public
** License version 2.1 as published by the Free Software Foundation and
** appearing in the file LICENSE.LGPL included in the packaging of this
** file. Please review the following information to ensure the GNU Lesser
** General Public License version 2.1 requirements will be met:
** http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** Other Usage
** Alternatively, this file may be used in accordance with the terms and
** conditions contained in a signed written agreement between you and Andrey Kartashov.
**
****************************************************************************/

#include "Arguments.hpp"
#include <QCoreApplication>

// Global static pointer used to ensure a single instance of the class.
Arguments* volatile Arguments::m_pArgumentsInstance = 0;
QSettings* Arguments::setup = 0;

QMutex Arguments::mutex;

/***************************************************************************************

****************************************************************************************/
QMap<QString,Arguments::_ArgDescr>& Arguments::getVarValStorage()
{
    static QMap<QString,_ArgDescr> m_VarValStorage;
    return m_VarValStorage;
}

/***************************************************************************************

****************************************************************************************/
Arguments& Arguments::Instance()
{
    if(!m_pArgumentsInstance)
    {
        mutex.lock();
        if(!m_pArgumentsInstance)
        {
            //Arguments* volatile temp =
            //      static_cast< Arguments* >(operator new (sizeof(Arguments)));
            m_pArgumentsInstance=new Arguments();
            setup=new QSettings(QSettings::IniFormat,QSettings::UserScope,QCoreApplication::organizationName());
            argsList();
        }
        mutex.unlock();
    }
    return *m_pArgumentsInstance;
}

/***************************************************************************************

****************************************************************************************/
Arguments::Arguments()
{
}

Arguments::~Arguments()
{
}


/***************************************************************************************

****************************************************************************************/
void Arguments::Init(QStringList l)
{
    int correct=0;
    foreach(const QString &key, getVarValStorage().keys())
    {
        int index=-1;
        if(getVarValStorage()[key]._cname!="")
        {
            if(getVarValStorage()[key]._type==QVariant::Bool)
            {
                index=l.indexOf(QRegExp("^\\-{1,2}"+getVarValStorage()[key]._cname+"$"),1);
            }
            else
            {
                index=l.indexOf(QRegExp("^\\-{1,2}"+getVarValStorage()[key]._cname+"=.+$"),1);
            }
            if(index>0)
            {
                QString command=l[index];

                if(command.left(2)=="--")
                {
                    command.remove(0,2);
                }
                else if(command.left(1)=="-")
                {
                    command.remove(0,1);
                }
                else
                {
                    usage();
                    throw "Incorrect parameter";
                }
                if(command=="help")
                {
                    usage();
                    exit(0);
                }
                /*Selecting correct command*/
                if(getVarValStorage()[key]._type==QVariant::Bool && command==getVarValStorage()[key]._cname)
                {
                    getVarValStorage()[key]._value=QVariant(true);
                    /*probably save to the INI*/
                    correct++;
                    continue;
                }
                if( getVarValStorage()[key]._type==QVariant::String &&
                        command.at(getVarValStorage()[key]._cname.size())==QChar('='))
                {
                    command.remove(0,getVarValStorage()[key]._cname.size()+1);
                    getVarValStorage()[key]._value=QVariant(command);
                    correct++;
                    continue;
                }
                if( getVarValStorage()[key]._type==QVariant::Int &&
                        command.at(getVarValStorage()[key]._cname.size())==QChar('='))
                {
                    command.remove(0,getVarValStorage()[key]._cname.size()+1);
                    getVarValStorage()[key]._value=QVariant(command).toInt();
                    correct++;
                    continue;
                }
                /*
                 *   by default treat all parameters as String
                 */
                if( command.at(getVarValStorage()[key]._cname.size())==QChar('='))
                {
                    command.remove(0,getVarValStorage()[key]._cname.size()+1);
                    getVarValStorage()[key]._value=QVariant(command);
                    correct++;
                    continue;
                }
                usage();
                throw "Incorrect parameter";
            }
            //index>0
        }

        if(!getVarValStorage()[key]._ininame.isEmpty() && setup->contains(getVarValStorage()[key]._ininame))
        {
            getVarValStorage()[key]._value=setup->value(getVarValStorage()[key]._ininame);
            continue;
        }

        if(!getVarValStorage()[key]._ininame.isEmpty() && getVarValStorage()[key]._stdin)
        {
            QTextStream stream(stdin);
            QString line;
#ifndef _WIN32
            termios oldt;
            tcgetattr(STDIN_FILENO, &oldt);
            termios newt = oldt;
            newt.c_lflag &= ~ECHO;
            tcsetattr(STDIN_FILENO, TCSANOW, &newt);
#else
            HANDLE hStdin = GetStdHandle(STD_INPUT_HANDLE);
            DWORD mode = 0;
            GetConsoleMode(hStdin, &mode);
            SetConsoleMode(hStdin, mode & (~ENABLE_ECHO_INPUT));
#endif
            cout<<"Please enter "<<getVarValStorage()[key]._descr.toStdString()<<":";

            line = stream.readLine();
            cout<<endl;
#ifndef _WIN32
            tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
#else
            SetConsoleMode(hStdin, mode & (ENABLE_ECHO_INPUT));
#endif
            QByteArray ba=qCompress(line.toLocal8Bit().data()).toBase64();
            getVarValStorage()[key]._value=ba;
            setup->setValue(getVarValStorage()[key]._ininame,ba);
            setup->sync();
            continue;
        }
        else
        {
            if(getVarValStorage()[key]._required==true)
            {
                usage();
                throw "Required parameter not present";
            }

        }
    }
    if(correct!=l.size()-1)
    {
        usage();
        for(int i=0;i<l.size();i++)
        {
            cout<<"DEBUG: "<<l[i].toStdString()<<endl;
        }
        throw "Incorrect command line";
    }
}

/***************************************************************************************

****************************************************************************************/
void Arguments::argsList(void)
{
    Arguments::addArg("in","in","inFileName",QVariant::String,"Input filename with MACS2 peak calling results, xls",QString(""));
    Arguments::addArg("out","out","outFileName",QVariant::String,"Base output file name",QString(""));
    Arguments::addArg("annotation","a","annotationFileName",QVariant::String,"Tab-separated annotation file",QString(""));
    Arguments::addArg("log","log","logFileName",QVariant::String,"log file name (default is ./logfile_def.log)",QString("./logfile_def.log"));
    Arguments::addArg("promoter","promoter","",QVariant::Int,"Promoter region around TSS in bp.",1000);
    Arguments::addArg("upstream","upstream","",QVariant::Int,"Upstream region before promoter in BP",20000);
    Arguments::addArg("sam_ignorechr","sam_ignorechr","SAM/IGNORECHR",QVariant::String,"The chromosome to be ignored",QString(""));//chrM
}

/***************************************************************************************

****************************************************************************************/
QFileInfo Arguments::fileInfo(const QString& str)
{
    return QFileInfo(gArgs().getArgs(str).toString());
}

/***************************************************************************************

****************************************************************************************/
QStringList  Arguments::split(const QString& str,const QChar& sep)
{
    return gArgs().getArgs(str).toString().split(sep);
}

/***************************************************************************************

****************************************************************************************/
void Arguments::usage(void)
{
   cout << "Usage:" <<endl;
   foreach(const QString &key,getVarValStorage().keys())
    {
        if(getVarValStorage()[key]._cname!="")
        {
            cout<<QString("\t--%1%2").arg(getVarValStorage()[key]._cname).
                  arg("",25-getVarValStorage()[key]._cname.length(),QChar(' ')).toStdString()<<"\t"<<getVarValStorage()[key]._descr.toStdString()<<endl;
        }
    }
}

/***************************************************************************************

****************************************************************************************/
void Arguments::addArg(QString key,QString _c/*command line argument*/, QString _i/*name in ini file*/,QVariant::Type _t, QString _d,QVariant _def,bool _r/*is argument required or not*/,bool _s/*should we read from stdin or not*/)
{
    if(!getVarValStorage().contains(key))
        getVarValStorage().insert(key,_ArgDescr(_c,_i,_t,_def,_d,_r,_s));
}

/***************************************************************************************

****************************************************************************************/
QVariant &Arguments::getArgs(QString key)
{
    mutex.lock();
    if(getVarValStorage().contains(key))
    {
        if(getVarValStorage()[key]._value.isValid() && !getVarValStorage()[key]._value.isNull())
        {
            QVariant &_v=getVarValStorage()[key]._value;
            mutex.unlock();
            return _v;
        }
        if(getVarValStorage()[key]._defValue.isValid() && !getVarValStorage()[key]._defValue.isNull())
        {
            QVariant &_v=getVarValStorage()[key]._defValue;
            mutex.unlock();
            return _v;
        }
    }

    mutex.unlock();
    throw "Incorrect key";
}
/***************************************************************************************

****************************************************************************************/

