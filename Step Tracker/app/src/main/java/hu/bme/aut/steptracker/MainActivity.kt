package hu.bme.aut.steptracker

import android.Manifest
import android.app.AlertDialog
import android.content.Context
import android.content.pm.PackageManager
import android.hardware.Sensor
import android.hardware.SensorEvent
import android.hardware.SensorEventListener
import android.hardware.SensorManager
import android.os.Bundle
import android.widget.Toast
import androidx.appcompat.app.AppCompatActivity
import androidx.appcompat.app.AppCompatDelegate
import androidx.core.app.ActivityCompat
import androidx.core.content.ContextCompat
import com.google.android.material.snackbar.Snackbar
import hu.bme.aut.steptracker.databinding.ActivityMainBinding
import hu.bme.aut.steptracker.model.Score
import hu.bme.aut.steptracker.sqlite.PersistentDataHelper
import java.util.*


class MainActivity : AppCompatActivity(), SensorEventListener {

    private lateinit var bindig: ActivityMainBinding
    private lateinit var sensorManager: SensorManager
    private var stepSensor: Sensor ?= null
    private var previousTotalSteps = 0f
    private lateinit var dataHelper: PersistentDataHelper
    private lateinit var today: String
    private val REQUEST_CODE = 102

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        bindig = ActivityMainBinding.inflate(layoutInflater)
        setContentView(bindig.root)
        dataHelper = PersistentDataHelper(this)
        sensorManager = getSystemService(SENSOR_SERVICE) as SensorManager
        today = getDate()
        readData()
        checkTheme()
        setupPermissions()
    }

    private fun setupPermissions() {
        if (ContextCompat.checkSelfPermission(this,
                Manifest.permission.ACTIVITY_RECOGNITION)
            != PackageManager.PERMISSION_GRANTED) {

            if (ActivityCompat.shouldShowRequestPermissionRationale(this,
                    Manifest.permission.ACTIVITY_RECOGNITION)) {
                Snackbar.make(bindig.root, "Permission needed for counting the steps!", Snackbar.LENGTH_LONG).show()
                makeRequest()
            } else {
                makeRequest()
            }
        } else {
            // Permission has already been granted
        }
    }

    private fun makeRequest() {
        ActivityCompat.requestPermissions(this,
            arrayOf(Manifest.permission.ACTIVITY_RECOGNITION),
            REQUEST_CODE)
    }

   override fun onRequestPermissionsResult(requestCode: Int, permissions: Array<String>,
   grantResults: IntArray) {
       super.onRequestPermissionsResult(requestCode, permissions, grantResults)
       when (requestCode) {
           REQUEST_CODE -> {
            if (grantResults.isNotEmpty() && grantResults[0] ==
            PackageManager.PERMISSION_GRANTED) {
               if ((ContextCompat.checkSelfPermission(this,
                       Manifest.permission.ACTIVITY_RECOGNITION)
                           != PackageManager.PERMISSION_GRANTED)) {
                   Snackbar.make(bindig.root, "Permission Granted", Snackbar.LENGTH_SHORT).show()
               }
            } else {
                AlertDialog.Builder(this)
                    .setTitle("Permission Needed")
                    .setMessage("Permission is needed to count your steps...\nChange it in your setting or reinstall the app")
                    .setPositiveButton("OK"
                    ) { _, _ ->
                        ActivityCompat.requestPermissions(
                            this@MainActivity,
                            arrayOf(Manifest.permission.READ_EXTERNAL_STORAGE),
                            REQUEST_CODE
                        )
                    }
                    .create().show()
            }
            return
         }
      }
   }

    private fun checkTheme() {
        val sharedPreference =  getSharedPreferences("Theme", Context.MODE_PRIVATE)
        val storedTheme: String? = sharedPreference?.getString("Theme","defaultName")
        if (storedTheme!=null && storedTheme =="Night"){
            AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_YES)
        }
    }

    override fun onResume() {
        super.onResume()
        stepSensor = sensorManager.getDefaultSensor(Sensor.TYPE_STEP_DETECTOR)

        if (stepSensor == null) {
            Toast.makeText(this, "No Step Counter Sensor Detected!", Toast.LENGTH_SHORT).show()
        }
        else {
            sensorManager.registerListener(this, stepSensor, SensorManager.SENSOR_DELAY_FASTEST)
        }
    }

    override fun onSensorChanged(event: SensorEvent?) {
        val newSteps = event!!.values[0]
        if (today == getDate()) {
            previousTotalSteps += newSteps
        }
        else {
            writeData(today)
            previousTotalSteps = newSteps
        }
        changeCounter()
        val sharedPreference =  getSharedPreferences("NumberOfSteps", Context.MODE_PRIVATE)
        val editor = sharedPreference.edit()
        editor.putFloat("steps",  previousTotalSteps)
        editor.apply()
        today = getDate()
    }

    override fun onAccuracyChanged(sensor: Sensor?, accuracy: Int) { }

    fun changeCounter () {
        bindig.stepValue.text = previousTotalSteps.toInt().toString()
    }

    override fun onStop() {
        super.onStop()
        writeData(today)
    }

    fun readData() {
        dataHelper.open()
        val score: Score? = dataHelper.searchScore(getDate())
        if (score != null){
            previousTotalSteps = score.steps
            changeCounter()
        }
        dataHelper.close()
    }

    fun writeData(date: String) {
        dataHelper.open()
        val score: Score? = dataHelper.searchScore(date)
        if (score == null)
            dataHelper.persistScore(Score(date, previousTotalSteps))
        else
            dataHelper.updateScore(date, previousTotalSteps)
        dataHelper.close()
    }

    fun getDate() : String {
        val c = Calendar.getInstance()
        return c.time.toString().dropLast(24) + c.time.toString().drop(29) + " " + c.get(Calendar.WEEK_OF_YEAR).toString()
    }
}