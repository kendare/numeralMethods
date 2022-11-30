// labwork1.cpp 

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>


//задание ТЕСТ
double func (double x, double y) // объявляю функцию задачи 1
{
    double f1 = 2 * x * exp(x) + y;
    return f1;
}
double funcT(double x, double y) // тестовой задачи
{
    double f1 = 0.16 * y;
    return f1;
}
std:: vector<double> RK_II_test(double x, double y, double h)
{
    std::vector<double> k{ 0, 0, 0 };
    k[1] = funcT(x, y);
    k[2] = funcT(x + h / 2, y + (h / 2) * k[1]);
    k[0] = y + h * k[2];
    return k;
}
std::vector<double> S2_test(double x, double y, double h)
{
    double y_i = y;
    double x_i = x;
    double h2 = h / 2;
    std::vector<double> vector{ 0, 0, 0, 0 }; // Первое значение S, второе y_i, третье y - y_i, четвертое  
    double v1 = RK_II_test(x, y, h)[0];
    for (int j = 0; j < 2; j++)
    {
        y_i = RK_II_test(x_i, y_i, h2)[0];
        x_i = x_i + h2;
    }
    vector[0] = abs((y_i - v1) / (pow(2, 2) - 1)); //s
    vector[1] = y_i;                                //y двойным шагом
    vector[2] = vector[0] * pow(2, 2);                 // LEE
    vector[3] = v1 - y_i;                      //численное - двойным шагом
    return vector;
}
void RK2_test()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector <double> s{ 0, 0, 0, 0 };
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ifstream fin;
    std::ofstream fout;
    const double eul = std::exp(1.0);
    double accuracy = 0.001;          // задаём точность
    double n = accuracy / pow(2, 3);
    double S = 0;
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    fout.open("test_RK2_output.txt");
    fin.open("test_input.txt");
    fin >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х
    double maxh = h;
    double minh = h;

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k2 = RK_II_test(xvector[i], yvector[i], h);
        y_1 = k2[0];
        s = S2_test(xvector[i], yvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        S = s[0];
        LEE = abs(s[2]);
        if (S < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S < n)
            {
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                c2++;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step
        double exactSol = pow(eul, 4 * xvector[i+1] / 25);
        /*std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k2[1] << "    k2 = "
            << k2[2] << "     LEE = " << LEE << "\n";*/
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << s[1]
            << "    \t yi - yi(1/2) = " << s[3]
            << "    \t LEE = " << LEE
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            <<"     \t y(xi) = " << exactSol
            <<"     \t |y(xi) - yi| = " << abs(exactSol - yvector[i+1])
            << "\n";
    }
    // std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = "<< c1 << "\nc2 = " << c2 << "\n\n";
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

std::vector<double> RK_III_test(double x, double y, double h)
{
    std::vector<double> k{ 0, 0, 0, 0, 0 };
    k[1] = funcT(x, y);
    k[2] = funcT(x + h / 2, y + h / 2 * k[1]);
    k[3] = funcT(x + h, y - h * k[1] + 2 * h * k[2]);
    k[0] = y + (h / 6) * (k[1] + 4 * k[2] + k[3]);
    return k;
}
std::vector <double> S3_test(double x, double y, double h)
{
    double y_i = y;
    double x_i = x;
    double h2 = h / 2;
    std::vector<double> vector{ 0, 0, 0, 0 }; // Первое значение S, второе y_i, третье y - y_i, четвертое  
    double v1 = RK_III_test(x, y, h)[0];
    for (int j = 0; j < 2; j++)
    {
        y_i = RK_III_test(x_i, y_i, h2)[0];
        x_i = x_i + h2;
    }
    vector[0] = abs((y_i - v1) / (pow(2, 3) - 1)); // s 
    vector[1] = y_i;                //половинным шагом
    vector[2] = vector[0] * pow(2, 3);//LEE
    vector[3] = v1 - y_i;
    return vector;
}
void RK3_test()
{

    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector <double> s{ 0, 0, 0, 0 };
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ofstream fout;
    std::ifstream fin;
    const double eul = std::exp(1.0);
    double accuracy = 0.0001;          // задаём точность
    double n = accuracy / pow(2, 4);
    double S = 0;
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    double maxh = h;
    double minh = 10;
    fout.open("test_RK3_output.txt");
    fin.open("test_input.txt");
    fin >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k3 = RK_III_test(xvector[i], yvector[i], h);
        y_1 = k3[0];
        s = S3_test(xvector[i], yvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        S = s[0];
        LEE = abs(s[2]);
        if (S < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S < n)
            {
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                c2++;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step
        double exactSol = pow(eul, 4 * xvector[i + 1] / 25);
        /*        std::cout << "n = " << i << "     y = "
                    << yvector[i + 1] << "    x = "
                    << xvector[i + 1] << "    h = " << h
                    << "   k1 = " << k3[1] << "    k2 = "
                    << k3[2] << "    k3 =" << k3[3] << "     LEE = " << LEE << "\n";*/                                  //тут вывод в консоль всех коэффициентов к и значений численных
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << s[1]
            << "    \t yi - yi(1/2) = " << s[3]
            << "    \t LEE = " << LEE
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "     \t y(xi) = " << exactSol
            << "     \t |y(xi) - yi| = " << abs(exactSol - yvector[i + 1])
            << "\n";
    }
    //std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();

}

std::vector < double > RK_IV_test(double x, double y, double h)
{
    std::vector<double> k4{ 0, 0, 0, 0, 0 };
    k4[1] = funcT(x, y);
    k4[2] = funcT(x + 0.5 * h, y + 0.5 * h * k4[1]);
    k4[3] = funcT(x + 0.5 * h, y + 0.5 * h * k4[2]);
    k4[4] = funcT(x + h, y + h * k4[3]);
    k4[0] = y + (h / 6) * (k4[1] + 2 * k4[2] + 2 * k4[3] + k4[4]);
    return k4;
}
std::vector <double> S4_test(double x, double y, double h)
{
    double y_i = y;
    double x_i = x;
    double h2 = h / 2; //считаем половинным шагом
    double v1 = 0;
    double v_1 = 0;
    std::vector<double> k1{ 0, 0, 0, 0, 0 };
    std::vector<double> k2{ 0, 0, 0, 0, 0 };
    std::vector<double> S{ 0, 0, 0 };
    for (int j = 0; j <= 1; j++)
    {
        k2 = RK_IV_test(x_i, y_i, h2); //k2 потому что двойной шаг
        x_i = x_i + h2;
        y_i = y_i + (h2 / 6) * (k2[1] + 2 * k2[2] + 2 * k2[3] + k2[4]);
    }
    S[0] = y_i;
    k1 = RK_IV_test(x, y, h);
    v1 = y + (h / 6) * (k1[1] + 2 * k1[2] + 2 * k1[3] + k1[4]);
    S[1] = (y_i - v1);
    S[2] = (S[1] / (pow(2, 4) - 1));
    return S;
}
void RK4_test()
{

    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector<double> S{ 0, 0, 0 };
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.00001 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счетчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ofstream fout;
    std::ifstream fin;
    const double eul = std::exp(1.0);
    double accuracy4 = 0.00001;          // задаём точность
    double n = accuracy4 / pow(2, 5);
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    double maxh = 0;
    double minh = 10;

    //здесь нужно получить все параметры h, x0, y0, b (правая граница)



    fout.open("test_RK4_output.txt");
    fin.open("test_input.txt");
    fin >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k4 = RK_IV_test(xvector[i], yvector[i], h);
        y_1 = k4[0];
        S = S4_test(xvector[i], yvector[i], h);
        S[1] = abs(S[1]);
        LEE = S[2] * pow(2, 4);
        if (S[2] < accuracy4)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S[2] < n)
            {
                c2++;
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на 
            }
        }
        else
        {
            c1++;
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            continue;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step

        double exactSol = pow(eul, 4 * xvector[i + 1] / 25);

        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << S[0]
            << "    \t yi - yi(1/2) = " << S[1]
            << "    \t LEE = " << LEE
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "    \t y(xi) = " << exactSol
            << "    \t |y(xi) - yi| = " << abs(exactSol - yvector[i + 1])
            << "\n";
    }
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

//Задание 1 2хе^x+y = y'
std::vector<double> RK_II(double x, double y, double h)
{
    std::vector<double> k{ 0, 0, 0 };
    k[1] = func(x, y);
    k[2] = func(x + h/2, y + ( h/2 )* k[1]);
    k[0] = y + h * k[2];
    return k;
}
std::vector<double> S2(double x, double y, double h)
{
    double y_i = y;
    double x_i = x;
    double h2 = h / 2;
    std::vector<double> vector{ 0, 0, 0, 0}; // Первое значение S, второе y_i, третье y - y_i, четвертое  
    double v1 = RK_II(x, y, h)[0];
    for (int j = 0; j < 2; j++)
    {
        y_i = RK_II(x_i, y_i, h2)[0];
        x_i = x_i + h2;
    }
    vector[0] = abs((y_i - v1) / (pow(2, 2) - 1)); //s
    vector[1] = y_i;                                //y двойным шагом
    vector[2] = vector[0]*pow(2,2);                 // LEE
    vector[3] = v1 - y_i;                      //численное - двойным шагом
    return vector;
}
void RK2()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector <double> s{ 0, 0, 0, 0 };
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ifstream fin;
    std::ofstream fout;
    double accuracy = 0.00000001;          // задаём точность
    double n = accuracy / pow(2, 3);
    double S = 0;
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    double maxh = h;
    double minh = 10;
    fout.open("RK2_output.txt");
    fin.open("input.txt");
    fin >> x0 >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k2 = RK_II(xvector[i], yvector[i], h);
        y_1 = k2[0];
        s = S2(xvector[i], yvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        S = s[0];
        LEE = abs(s[2]);
        if (S < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S < n)
            {
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                c2++;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step
        /*std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k2[1] << "    k2 = "
            << k2[2] << "     LEE = " << LEE << "\n";*/
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << s[1]
            << "    \t yi - yi(1/2) = " << s[3]
            << "    \t LEE = " << LEE
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "\n";
    }
   // std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = "<< c1 << "\nc2 = " << c2 << "\n\n";
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

std::vector<double> RK_III(double x, double y, double h)
{
    std::vector<double> k{ 0, 0, 0, 0, 0 };
    k[1] = func(x, y);
    k[2] = func(x + h / 2, y + h / 2 * k[1]);
    k[3] = func(x + h, y - h * k[1] + 2 * h * k[2]);
    k[0] = y + (h / 6) * (k[1] + 4 * k[2] + k[3]);
    return k;
}
std::vector<double> S3(double x, double y, double h)
{
    double y_i = y;
    double x_i = x;
    double h2 = h / 2;
    std::vector<double> vector{ 0, 0, 0, 0 }; // Первое значение S, второе y_i, третье y - y_i, четвертое  
    double v1 = RK_III(x, y, h)[0];
    for (int j = 0; j < 2; j++)
    {
        y_i = RK_III(x_i, y_i, h2)[0];
        x_i = x_i + h2;
    }
    vector[0] = abs((y_i - v1) / (pow(2, 3) - 1)); // s 
    vector[1] = y_i;                //половинным шагом
    vector[2] = vector[0]*pow(2, 3);//LEE
    vector[3] = v1 - y_i;
    return vector;
}
void RK3()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector <double> s{0, 0, 0, 0};
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ofstream fout;
    std::ifstream fin;
    double accuracy = 0.00000001;          // задаём точность
    double n = accuracy / pow(2, 4);
    double S = 0;
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    double maxh = h;
    double minh = 10;
    fout.open("RK3_output.txt");
    fin.open("input.txt");
    fin >> x0 >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k3 = RK_III(xvector[i], yvector[i], h);
        y_1 = k3[0];
        s = S3(xvector[i], yvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        S = s[0];
        LEE = abs(s[2]);
        if (S < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S < n)
            {
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                c2++;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step
/*        std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k3[1] << "    k2 = "
            << k3[2] << "    k3 =" << k3[3] << "     LEE = " << LEE << "\n";*/                                  //тут вывод в консоль всех коэффициентов к и значений численных
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << s[1]
            << "    \t yi - yi(1/2) = " << s[3]
            << "    \t LEE = " << LEE
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "\n";
    }
    //std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

std::vector<double> RK_IV(double x, double y, double h) //Считаю коэффициенты для 
{
    std::vector<double> k4{0, 0, 0, 0, 0};
    k4[1] = func(x, y);
    k4[2] = func(x + 0.5 * h, y + 0.5 * h * k4[1]);
    k4[3] = func(x + 0.5 * h, y + 0.5 * h *  k4[2]);
    k4[4] = func(x + h, y + h * k4[3]);
    k4[0] = y +(h / 6) * (k4[1] + 2 * k4[2] + 2 * k4[3] + k4[4]);
    return k4;
}
std::vector<double> S4(double xi, double yi, double h)
{
    double y_i = yi;
    double x_i = xi;
    double h2 = h / 2; //считаем половинным шагом
    double v1 = 0;
    double v_1 = 0;
    std::vector<double> k1{ 0, 0, 0, 0, 0 };
    std::vector<double> k2{ 0, 0, 0, 0, 0 };
    std::vector<double> S{ 0, 0, 0};
    for (int j = 0; j <= 1; j++)
    {
        k2 = RK_IV(x_i, y_i, h2); //k2 потому что двойной шаг
        x_i = x_i + h2;
        y_i = y_i + (h2 / 6) * (k2[1] + 2 * k2[2] + 2 * k2[3] + k2[4]);
    }
    S[0] = y_i;
    k1 = RK_IV(xi, yi, h);
    v1 = yi + (h / 6) * (k1[1] + 2 * k1[2] + 2 * k1[3] + k1[4]);
    S[1] = (y_i - v1);
    S[2] = (S[1] / (pow(2, 4) - 1));
    return S;
}
void RK4()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector<double> S{0, 0, 0};
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double h{ 0.00001 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счетчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ofstream fout;
    std::ifstream fin;
    double accuracy4 = 0.000000001;          // задаём точность
    double n = accuracy4 / pow(2, 5);
    double x_1 = 0;
    double y_1 = 0;
    double LEE = 0;
    double maxLEE = 0;
    double maxh = 0;
    double minh = 10;

    //здесь нужно получить все параметры h, x0, y0, b (правая граница)



    fout.open("RK4_output.txt");
    fin.open("input.txt");
    fin >> x0 >> y0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        std::vector<double> k4 = RK_IV(xvector[i], yvector[i], h);
        y_1 = k4[0];
        S = S4(xvector[i], yvector[i], h);
        S[1] = abs(S[1]);
        LEE = S[2] * pow(2, 4);
        if (S[2] < accuracy4)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            yvector.push_back(y_1);         //Здесь скорее всего узкое место, нужно посмотреть как добавить новый элемент
            xvector.push_back(x_1);
            if (S[2] < n)
            {
                c2++;
                h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на 
            }
        }
        else
        {
            c1++;
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            continue;
        }
        if (LEE > maxLEE) maxLEE = LEE;     // max Local accuracy
        if (h > maxh) maxh = h;             // max step
        if (h < minh) minh = h;             // min step

        fout <<"\t i = " << i 
             <<"    \t yi = " << yvector[i + 1]
             <<"    \t xi = " << xvector[i + 1]
             <<"    \t h = " << h 
             <<"    \t yi(1/2) = "<< S[0]
             <<"    \t yi - yi(1/2) = "<< S[1]
             <<"    \t LEE = " << LEE
             <<"    \t c1 = " << c1
             <<"    \t c2 = " << c2
             << "\n";
    }
    fout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}


//Задание 2 y' = x^2 - z, z' = y+x, y(0) = 1, z(0) = 1, b = 1;

double funcY(double x, double y, double z) // функция y
{
    double y_ = pow(x, 2) - z;
    return y_;
}

double funcZ(double x, double y, double z) //функция z
{
    double z_ = y + x;
    return z_;
}

std::vector<double> RK_II_sys(double x, double y, double z, double h)
{
    std::vector<double> result{ 0,0 };
    std::vector<double> k{0, 0, 0, 0};
    k[0] = funcY(x, y, z);
    k[1] = funcZ(x, y, z);
    k[2] = funcY(x + h / 2, y + (h / 2) * k[0], z + (h / 2) * k[1]);
    k[3] = funcZ(x + h / 2, y + (h / 2) * k[0], z + (h / 2) * k[1]);
    result[0] = y + h * k[2];        //реузльтирующая у
    result[1] = z + h * k[3];        //результирующая z
    return result;
}
std::vector<double> S2_sys(double x, double y, double z, double h)
{
    std::vector<double> S{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<double> result{ 0, 0 };
    double yi = y;
    double zi = z;
    double xi = x;
    double h2 = h / 2;
    for (int i = 0; i < 2; i++)
    {
        result = RK_II_sys(xi, yi, zi, h2);
        xi = xi + h2;
        yi = result[0];
        zi = result[1];
    }
    std::vector<double> preResult{ 0,0 };
    preResult = RK_II_sys(x, y, z, h);
    S[0] = preResult[0];// y
    S[1] = preResult[1];// z
    S[3] = yi;          // y 1/2
    S[4] = zi;          // z 1/2
    S[5] = S[0] - S[3]; // y - y1/2
    S[6] = S[1] - S[4]; // z - z1/2
    S[7] = abs(S[5]) / (pow(2, 2) - 1); // Sy
    S[8] = abs(S[7]) / (pow(2, 2) - 1); // Sz
    S[9] = S[7] * pow(2, 2);            //LEEy
    S[10] = S[8] * pow(2, 2);           //LEEZ
    return S;
}
void RK2_sys()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector<double> zvector;
    std::vector<double> S(11);
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double z0{ 1 };
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ifstream fin;
    std::ofstream fout;
    double accuracy = 0.001;          // задаём точность
    double n = accuracy / pow(2, 3);
    double x_1 = 0;
    double y_1 = 0;
    double z_1 = 0;
    fout.open("sys_RK2_output.txt");
    fin.open("sys_input.txt");
    fin >> x0 >> y0 >> z0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х
    zvector.push_back(z0);  // Добавили z0 к значениям вектора z

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        S = S2_sys(xvector[i], yvector[i], zvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        if (S[7] < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            if (S[8] < accuracy)
            {
                yvector.push_back(S[0]);
                zvector.push_back(S[1]);
                xvector.push_back(x_1);
                if (S[7] < n)
                    if (S[8] < n)
                    {
                        h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                        c2++;
                    }
            }
            else
            {
                h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
                i = i - 1;
                c1++;
                continue;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
            continue;
        }
        /*std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k2[1] << "    k2 = "
            << k2[2] << "     LEE = " << LEE << "\n";*/
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t zi = " << zvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << S[3]
            << "    \t zi(1/2) = " << S[4]
            << "    \t yi - yi(1/2) = " << S[5]
            << "    \t zi - zi(1/2) = " << S[6]
            << "    \t LEEy = " << S[9]
            << "    \t LEEy = " << S[10]
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "\n";
    }
    // std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = "<< c1 << "\nc2 = " << c2 << "\n\n";
    fout << "\n\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

std::vector<double> RK_III_sys(double x, double y, double z, double h)
{
    std::vector<double> result{ 0, 0 };
    std::vector<double> k(6);
    k[0] = funcY(x, y, z);
    k[1] = funcZ(x, y, z);
    k[2] = funcY(x + h / 2, y + h / 2 * k[0], z + h/2 * k[1]);
    k[3] = funcZ(x + h / 2, y + h / 2 * k[0], z + h / 2 * k[1]);
    k[4] = funcY(x + h, y - h * k[0] + 2 * h * k[2], z - h * k[1] + 2 * h * k[3]);
    k[5] = funcZ(x + h, y - h * k[0] + 2 * h * k[2], z - h * k[1] + 2 * h * k[3]);
    result[0] = y + (h / 6) * (k[0] + 4 * k[2] + k[4]);
    result[1] = z + (h / 6) * (k[1] + 4 * k[3] + k[5]);
    return result;
}
std::vector<double> S3_sys(double x, double y, double z, double h)
{
    std::vector<double> S(11);
    std::vector<double> result{ 0, 0 };
    double yi = y;
    double zi = z;
    double xi = x;
    double h2 = h / 2;
    for (int i = 0; i < 2; i++)
    {
        result = RK_III_sys(xi, yi, zi, h2);
        xi = xi + h2;
        yi = result[0];
        zi = result[1];
    }
    std::vector<double> preResult{ 0,0 };
    preResult = RK_III_sys(x, y, z, h);
    S[0] = preResult[0];// y
    S[1] = preResult[1];// z
    S[3] = yi;          // y 1/2
    S[4] = zi;          // z 1/2
    S[5] = S[0] - S[3]; // y - y1/2
    S[6] = S[1] - S[4]; // z - z1/2
    S[7] = abs(S[5]) / (pow(2, 2) - 1); // Sy
    S[8] = abs(S[7]) / (pow(2, 2) - 1); // Sz
    S[9] = S[7] * pow(2, 2);            //LEEy
    S[10] = S[8] * pow(2, 2);           //LEEZ
    return S;
}
void RK3_sys()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector<double> zvector;
    std::vector<double> S(11);
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double z0{ 1 };
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ifstream fin;
    std::ofstream fout;
    double accuracy = 0.0001;          // задаём точность
    double n = accuracy / pow(2, 3);
    double x_1 = 0;
    double y_1 = 0;
    double z_1 = 0;
    fout.open("sys_RK3_output.txt");
    fin.open("sys_input.txt");
    fin >> x0 >> y0 >> z0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х
    zvector.push_back(z0);  // Добавили z0 к значениям вектора z

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        S = S3_sys(xvector[i], yvector[i], zvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        if (S[7] < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            if (S[8] < accuracy)
            {
                yvector.push_back(S[0]);
                zvector.push_back(S[1]);
                xvector.push_back(x_1);
                if (S[7] < n)
                    if (S[8] < n)
                    {
                        h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                        c2++;
                    }
            }
            else
            {
                h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
                i = i - 1;
                c1++;
                continue;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
            continue;
        }
        /*std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k2[1] << "    k2 = "
            << k2[2] << "     LEE = " << LEE << "\n";*/
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t zi = " << zvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << S[3]
            << "    \t zi(1/2) = " << S[4]
            << "    \t yi - yi(1/2) = " << S[5]
            << "    \t zi - zi(1/2) = " << S[6]
            << "    \t LEEy = " << S[9]
            << "    \t LEEy = " << S[10]
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "\n";
    }
    // std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = "<< c1 << "\nc2 = " << c2 << "\n\n";
    fout << "\n\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

std::vector<double> RK_IV_sys(double x, double y, double z, double h)
{
    std::vector<double> result(2);
    std::vector<double> k(8);
    k[0] = funcY(x, y, z);
    k[1] = funcZ(x, y, z);
    k[2] = funcY(x + 0.5 * h, y + 0.5 * h * k[0], z + h/2 * k[1]);
    k[3] = funcZ(x + 0.5 * h, y + 0.5 * h * k[0], z + h / 2 * k[1]);
    k[4] = funcY(x + 0.5 * h, y + 0.5 * h * k[2], z + h / 2 * k[3]);
    k[5] = funcZ(x + 0.5 * h, y + 0.5 * h * k[2], z + h / 2 * k[3]);
    k[6] = funcY(x + h, y + h * k[4], z + h * k[5]);
    k[7] = funcZ(x + h, y + h * k[4], z + h * k[5]);

    result[0] = y + (h / 6) * (k[0] + 2 * k[2] + 2 * k[4] + k[6]);
    result[1] = z + (h / 6) * (k[1] + 2 * k[3] + 2 * k[5] + k[7]);

    return result;
}
std::vector<double> S4_sys(double x, double y, double z, double h)
{
    std::vector<double> S(11);
    std::vector<double> result{ 0, 0 };
    double yi = y;
    double zi = z;
    double xi = x;
    double h2 = h / 2;
    for (int i = 0; i < 2; i++)
    {
        result = RK_IV_sys(xi, yi, zi, h2);
        xi = xi + h2;
        yi = result[0];
        zi = result[1];
    }
    std::vector<double> preResult{ 0,0 };
    preResult = RK_IV_sys(x, y, z, h);
    S[0] = preResult[0];// y
    S[1] = preResult[1];// z
    S[3] = yi;          // y 1/2
    S[4] = zi;          // z 1/2
    S[5] = S[0] - S[3]; // y - y1/2
    S[6] = S[1] - S[4]; // z - z1/2
    S[7] = abs(S[5]) / (pow(2, 2) - 1); // Sy
    S[8] = abs(S[7]) / (pow(2, 2) - 1); // Sz
    S[9] = S[7] * pow(2, 2);            //LEEy
    S[10] = S[8] * pow(2, 2);           //LEEZ
    return S;
}
void RK4_sys()
{
    std::vector<double> yvector;
    std::vector<double> xvector;
    std::vector<double> zvector;
    std::vector<double> S(11);
    double x0{ 0 };     // объявляем х0
    double y0{ 1 };     // объявляем у0
    double z0{ 1 };
    double h{ 0.1 };    // объявляем шаг h
    int b{ 1 };         // правая граница
    int c1{ 0 };        // счётчик делений шага
    int c2{ 0 };        // счетчик удвоений шага
    std::ifstream fin;
    std::ofstream fout;
    double accuracy = 0.000001;          // задаём точность
    double n = accuracy / pow(2, 3);
    double x_1 = 0;
    double y_1 = 0;
    double z_1 = 0;
    fout.open("sys_RK4_output.txt");
    fin.open("sys_input.txt");
    fin >> x0 >> y0 >> z0 >> h >> b;
    yvector.push_back(y0);  //Добавили y0 к значениям вектора у
    xvector.push_back(x0);  //Добавили х0 к значениям вектора х
    zvector.push_back(z0);  // Добавили z0 к значениям вектора z

    for (int i = 0; xvector[i] < b; i++) //Считаем методом 4 порядка.
    {
        x_1 = xvector[i] + h;
        if (x_1 > b)
        {
            h = b - xvector[i];
            x_1 = b;
        }
        S = S4_sys(xvector[i], yvector[i], zvector[i], h);  // 0 - S, 1 - LEE, 2 - y половинным шагом 
        if (S[7] < accuracy)                  //если значение удовлетворяет условиям, тогда принимаем его
        {
            if (S[8] < accuracy)
            {
                yvector.push_back(S[0]);
                zvector.push_back(S[1]);
                xvector.push_back(x_1);
                if (S[7] < n)
                    if (S[8] < n)
                    {
                        h = h * 2;                  // проверяем не слишком ли маленький шаг. Если слишком маленький, то умножаем на два
                        c2++;
                    }
            }
            else
            {
                h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
                i = i - 1;
                c1++;
                continue;
            }
        }
        else
        {
            h = h / 2;                      //Если шаг слишком большой, то делим его пополам и повторяем цикл заново
            i = i - 1;
            c1++;
            continue;
        }
        /*std::cout << "n = " << i << "     y = "
            << yvector[i + 1] << "    x = "
            << xvector[i + 1] << "    h = " << h
            << "   k1 = " << k2[1] << "    k2 = "
            << k2[2] << "     LEE = " << LEE << "\n";*/
        fout << "\t i = " << i
            << "    \t yi = " << yvector[i + 1]
            << "    \t zi = " << zvector[i + 1]
            << "    \t xi = " << xvector[i + 1]
            << "    \t h = " << h
            << "    \t yi(1/2) = " << S[3]
            << "    \t zi(1/2) = " << S[4]
            << "    \t yi - yi(1/2) = " << S[5]
            << "    \t zi - zi(1/2) = " << S[6]
            << "    \t LEEy = " << S[9]
            << "    \t LEEy = " << S[10]
            << "    \t c1 = " << c1
            << "    \t c2 = " << c2
            << "\n";
    }
    // std::cout << "max LEE = " << maxLEE << "\nmax h = " << maxh << "\nmin h = " << minh << "\nc1 = "<< c1 << "\nc2 = " << c2 << "\n\n";
    fout << "\n\nc1 = " << c1 << "\nc2 = " << c2 << "\n\n";
    fout.close();
}

int main()
{
    setlocale(LC_ALL, "Russian");
    int z;
    int c;
    int a;
    int k;
    do
    {
        std::cout << "Выберите какую задачу вы хотите решить\n"
            << " 1 - Тестовая задача\n"
            << " 2 - Задача 1\n"
            << " 3 - Задача 2\n"
            << " 0 - Закрыть программу"
            << "\n\n\n";
        std::cin >> c;
        std::cout << "\n\n\n";
        switch (c)
        {
        case 1:
            do
            {
                std::cout << "Выберите каким методом вы хотите найти знаечния ТЕСТОВОЙ задачи\n";
                std::cout << "1  -  прогнать методом Рунге-Кутты 2 порядка\n";
                std::cout << "2  -  прогнать методом Рунге-Кутты 3 порядка\n";
                std::cout << "3  -  прогнать методом Рунге-Кутты 4 порядка\n";
                std::cout << "0  -  выйти из этого раздела меню\n\n" << std::endl;
                std::cout << "\n\n\n";
                std::cin >> a;
                std::cout << "\n\n\n";
                switch (a)
                {
                case 1:
                    RK2_test();
                    break;
                case 2:
                    RK3_test();
                    break;
                case 3:
                    RK4_test();
                    break;
                case 0:
                    a = 0;
                    break;
                default:
                    std::cout << "Введено некорректное значение\n" << std::endl;
                    break;
                }
            } while (a != 0);
            break;
        case 2:
            do
            {
                std::cout << "Выберите каким методом вы хотите найти знаечния задачи 1\n";
                std::cout << "1  -  прогнать методом Рунге-Кутты 2 порядка\n";
                std::cout << "2  -  прогнать методом Рунге-Кутты 3 порядка\n";
                std::cout << "3  -  прогнать методом Рунге-Кутты 4 порядка\n";
                std::cout << "0  -  выйти из этого раздела меню\n\n" << std::endl;
                std::cout << "\n\n\n";
                std::cin >> z;
                std::cout << "\n\n\n";
                switch (z)
                {
                case 1:
                {
                    RK2();
                    break;
                }
                case 2:
                {
                    RK3();
                    break;
                }
                case 3:
                {
                    RK4();
                    break;
                }
                case 0:
                    break;
                default:
                {
                    std::cout << "Введено некорректное значение\n" << std::endl;
                    break;
                }
                }
            } while (z != 0);
            break;
        case 3:
            do {
                std::cout << "Выберите каким методом вы хотите найти знаечния задачи 2\n";
                std::cout << "1  -  прогнать методом Рунге-Кутты 2 порядка\n";
                std::cout << "2  -  прогнать методом Рунге-Кутты 3 порядка\n";
                std::cout << "3  -  прогнать методом Рунге-Кутты 4 порядка\n";
                std::cout << "0  -  выйти из этого раздела меню\n\n" << std::endl;
                std::cout << "\n\n\n";
                std::cin >> k;
                std::cout << "\n\n\n";
                switch (k)
                {
                case 0:
                    break;
                case 1:
                    RK2_sys();
                    break;
                case 2:
                    RK3_sys();
                    break;
                case 3:
                    RK4_sys();
                    break;
                default:
                    std::cout << "Введите корректное значение";
                }
            } while (k != 0);
            break;
        case 0:
            break;
        default:
            std::cout << "Введите корректное значение";
            break;
        }
    } while (c != 0);
}