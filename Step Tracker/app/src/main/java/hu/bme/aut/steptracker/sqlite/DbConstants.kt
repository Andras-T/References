package hu.bme.aut.steptracker.sqlite

import android.database.sqlite.SQLiteDatabase
import android.util.Log

object DbConstants {
    const val DATABASE_NAME = "stepcounter.db"
    const val DATABASE_VERSION = 1

    object Score {
        const val DATABASE_TABLE = "steps"

        enum class Columns {
            ID, NUMBEROFSTEPS, DATE
        }

        private val DATABASE_CREATE = """create table if not exists $DATABASE_TABLE (
            ${Columns.ID.name} integer primary key autoincrement,
            ${Columns.NUMBEROFSTEPS.name} real,
            ${Columns.DATE} real not null
            );"""

        private const val DATABASE_DROP = "drop table if exists $DATABASE_TABLE;"

        fun onCreate(database: SQLiteDatabase) {
            database.execSQL(DATABASE_CREATE)
        }

        fun onUpgrade(database: SQLiteDatabase, oldVersion: Int, newVersion: Int) {
            Log.w(
                Score::class.java.name,
                "Upgrading from version $oldVersion to $newVersion"
            )
            database.execSQL(DATABASE_DROP)
            onCreate(database)
        }
    }
}