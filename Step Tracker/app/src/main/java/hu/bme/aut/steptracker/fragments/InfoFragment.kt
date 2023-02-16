package hu.bme.aut.steptracker.fragments

import android.content.Context
import android.os.Bundle
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.fragment.app.Fragment
import androidx.fragment.app.FragmentActivity
import hu.bme.aut.steptracker.databinding.FragmentInfoBinding
import hu.bme.aut.steptracker.model.Score
import hu.bme.aut.steptracker.sqlite.PersistentDataHelper


class InfoFragment : Fragment() {

    private lateinit var binding: FragmentInfoBinding
    private lateinit var dataHelper: PersistentDataHelper

    override fun onCreateView(inflater: LayoutInflater, container: ViewGroup?, savedInstanceState: Bundle?): View {
        super.onCreateView(inflater, container, savedInstanceState)
        binding = FragmentInfoBinding.inflate(inflater, container, false)
        dataHelper = activity?.let { PersistentDataHelper(it) }!!
        readData()
        return binding.root
    }

    private fun readData() {
        dataHelper.open()
        // Total steps
        var steps = 0
        val scores = dataHelper.getScore()
        var max = 0f
        for (score: Score in scores){
            steps += score.steps.toInt()
            if (score.steps > max)
                max = score.steps
        }
        binding.totalSteps.text = steps.toString()
        // Calories burned
        val act: FragmentActivity? = this.activity
        if (act != null) {
            val sharedPreferenceWeight =  act.getSharedPreferences("Weight", Context.MODE_PRIVATE)
            val weight = sharedPreferenceWeight.getString("Weight","75")
            val sharedPreferenceHeight =  act.getSharedPreferences("Height", Context.MODE_PRIVATE)
            val height = sharedPreferenceHeight.getString("Height","165")
            if (weight != null && height != null){
                var kcal = (weight.toInt() * height.toInt() * 0.0003154 * 0.01 * steps).toString()
                if(kcal.length > 6){
                    kcal = kcal.dropLast(kcal.length-6)
                }
                binding.caloriesBurned.text = kcal
            }
        }
        // Highest score
        binding.highestScore.text = max.toInt().toString()
        dataHelper.close()
    }
}