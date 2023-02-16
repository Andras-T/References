package hu.bme.aut.steptracker.sqlite

import android.content.Context
import android.database.Cursor
import android.database.sqlite.SQLiteDatabase
import android.database.sqlite.SQLiteOpenHelper
import java.util.*

//TimeFormat: Sat Nov 19 2022 47
class DbHelper(context: Context) :
    SQLiteOpenHelper(context, DbConstants.DATABASE_NAME, null, DbConstants.DATABASE_VERSION) {

    override fun onCreate(sqLiteDatabase: SQLiteDatabase) {
        DbConstants.Score.onCreate(sqLiteDatabase)
    }

    override fun onUpgrade(
        sqLiteDatabase: SQLiteDatabase,
        oldVersion: Int,
        newVersion: Int
    ) {
        DbConstants.Score.onUpgrade(sqLiteDatabase, oldVersion, newVersion)
    }

    fun getScore(date: String): Cursor {
        val db = this.readableDatabase
        return db.rawQuery("SELECT * FROM steps WHERE steps.DATE='$date' ", null)
    }

    fun updateScore(date: String, score: Float) {
        val db = this.writableDatabase
        val strSQL = "UPDATE steps SET NUMBEROFSTEPS = $score WHERE DATE ='$date'"
        db.execSQL(strSQL)
    }

    fun getWeeklyScore(week: String) : Cursor{
        val db = this.writableDatabase
        val year = Calendar.getInstance().get(Calendar.YEAR).toString()
        return db.rawQuery("SELECT * FROM steps WHERE steps.DATE LIKE '%$year $week'", null)
    }
}