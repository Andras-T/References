#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>

template <typename T> class vectorClass
{
private:

    //tarolo
    T* arr;
    // tarola kapacitasa
    size_t kapacitas;
    // aktualis meret
    size_t aktualis;

public:

    /*!
     *  Konstruktor.
     */
    vectorClass()
    {
        arr = 0;
        kapacitas = 0;
        aktualis = 0;
    }

    /*!
     *  Destruktor.
     */
    ~vectorClass() 
    {
        if (arr != NULL)
        {
            delete[] arr;
        }   
    }

    /*!
     *  Operator= vectorClass
     *
     *      @param masik vectorClass
     *
     *      @return vectorClass
     */
    vectorClass& operator=(const vectorClass& other)
    {
        if (kapacitas != 0)
        {
            delete[] arr;
        }
        aktualis = other.aktualis;
        kapacitas = other.kapacitas;
        if (other.arr!= 0)
        {
        arr = new T[kapacitas];
        for (size_t i = 0; i < aktualis; i++)
        {
            arr[i] = other.arr[i];
        }
        }
        return *this;
    }

    /*!
     *  uj elem hozzaadasa
     *
     *      @param adat 
     */ 
    void push(T data)
    {
        if (kapacitas == 0)
        {
            arr = new T[1];
            kapacitas = 1;
        }
        // tarolo bavitese 2*
        if (aktualis == kapacitas) 
        {
            T* temp = new T[2 * kapacitas];

            // regi elemek atmasolasa
            for (size_t i  = 0; i < kapacitas; i++) 
            {
                temp[i] = arr[i];
            }

            // regi tomb torlese
            delete[] arr;
            kapacitas *= 2;
            arr = temp;
        }
        
        // adatok visszahelyezese
        arr[aktualis] = data;
        aktualis++;
    }

    /*!
     *  indexelo operator
     *
     *      @param index 
     *
     *      @return Tarolo i. eleme
     */
    T& operator[](const size_t index) const
    {
        if (index >= 0 && index < aktualis)
            return arr[index];
        else
            throw std::invalid_argument("Rossz index");
    }

    /*!
     *  utolso elem eltavolitasa
     */
    void pop() 
    {
        if (aktualis > 0)
            aktualis--;
        else
            throw std::invalid_argument("Nincs tobb");
    }

    /*!
     *  beolvasott elemek szamat adja meg
     *
     *      @return beolvasott elemek
     */
    size_t size() const { return aktualis; }

    /*!
     *  szerepel-e a taroloban
     *
     *      @param adat
     *
     *      @return szerepel-e
     */
    bool init(T data) const
    {
        for (size_t i = 0; i < aktualis; i++)
        {
            if (data == arr[i])
                return true;
        }
        return false;
    }
};

//
/*! 
 *      inserter
 *      @param [in,out] os  
 *      @param [in]     vectorClass 
 *
 *      @return os
 */
template<class K> inline std::ostream& operator<<(std::ostream& os, const vectorClass<K>& val) {
    for (size_t i = 0; i < val.size(); i++)
    {
        os << val[i];
    }

    return os;
}

#endif //VECTOR_H
