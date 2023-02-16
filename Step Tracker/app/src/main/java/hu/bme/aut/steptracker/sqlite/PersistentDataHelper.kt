package hu.bme.aut.steptracker.sqlite

import android.content.ContentValues
import android.content.Context
import android.database.Cursor
import android.database.sqlite.SQLiteDatabase
import android.database.sqlite.SQLiteException
import hu.bme.aut.steptracker.model.Score

class PersistentDataHelper(context: Context) {
    private var database: SQLiteDatabase? = null
    private val dbHelper: DbHelper = DbHelper(context)

    private val scoreColumns = arrayOf(
        DbConstants.Score.Columns.ID.name,
        DbConstants.Score.Columns.NUMBEROFSTEPS.name,
        DbConstants.Score.Columns.DATE.name
    )


    @Throws(SQLiteException::class)
    fun open() {
        database = dbHelper.writableDatabase
    }

    fun close() {
        dbHelper.close()
    }

    fun searchScore(date: String) : Score?{
        val score: Score
        val cursor: Cursor = dbHelper.getScore(date)
        if (cursor.moveToFirst()) {
            score = cursorToScore(cursor)
        } 
        else {
            return null
        }
        return score
    }

    fun searchWeeklyScore(week: String) :MutableList<Score> {
        val scores: MutableList<Score> = ArrayList()
        val cursor: Cursor = dbHelper.getWeeklyScore(week)
        cursor.moveToFirst()
        while (!cursor.isAfterLast) {
            val score: Score = cursorToScore(cursor)
            scores.add(score)
            cursor.moveToNext()
        }
        cursor.close()
        return scores
    }

    fun updateScore(date: String, score: Float) {
        dbHelper.updateScore(date, score)
    }

    fun persistScore(score: Score) {
            val values = ContentValues()
            values.put(DbConstants.Score.Columns.NUMBEROFSTEPS.name, score.steps)
            values.put(DbConstants.Score.Columns.DATE.name, score.date)
            database!!.insert(DbConstants.Score.DATABASE_TABLE, null, values)
    }

    fun getScore(): MutableList<Score> {
        val scores: MutableList<Score> = ArrayList()
        val cursor: Cursor =
            database!!.query(DbConstants.Score.DATABASE_TABLE, scoreColumns, null, null, null, null, null)
        cursor.moveToFirst()
        while (!cursor.isAfterLast) {
            val score: Score = cursorToScore(cursor)
            scores.add(score)
            cursor.moveToNext()
        }
        cursor.close()
        return scores
    }

    private fun cursorToScore(cursor: Cursor): Score {
        return Score(
            cursor.getString(DbConstants.Score.Columns.DATE.ordinal),
            cursor.getFloat(DbConstants.Score.Columns.NUMBEROFSTEPS.ordinal)
        )
    }

}