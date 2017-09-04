#include <QCoreApplication>
#include <qDebug>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    qDebug() <<  "hello world";

    return a.exec();
}
