// stepik.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "md.h"

int main()
{
    test();

    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
    omp_set_num_threads(8);

    ssystem sample(force_LD, potential_LD);
    cout << "Okay, lets go!\n" << pow(2., 1. / 6.) << endl;
    int number_of_steps = 80000;                            // число МД шагов
    sample.initialize();                                  // инициализируем систему
    int time_for_kinE_count = 0.1 * number_of_steps;
    int counter_for_kinE = 0;
    double start = omp_get_wtime();
    for (int i = 0; i < number_of_steps; i++)
    {
        cout << "calc_forces\n";
        sample.calc_forces();
        cout << "integrate\n";
        sample.integrate();
        if (counter_for_kinE == time_for_kinE_count)
        {
            cout << "MD step #" << i << " time: " << (omp_get_wtime() - start) << endl;
            sample.measure_all();                         // вывод промежуточных итогов
            counter_for_kinE = 0;
            start = omp_get_wtime();
        }
        else
            counter_for_kinE++;

    }
    sample.output_data();

    cout << "The end!\n";
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
