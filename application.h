#ifndef APPLICATION_H
#define APPLICATION_H

#include <QApplication>
#include <QCoreApplication>
#include <QMutex>
#define App (static_cast<application*>(QCoreApplication::instance()))
class application : public QApplication
{
    Q_OBJECT
public:
    application(int& argc, char **argv);
    ~application();
signals:


public:
};

#endif // APPLICATION_H
