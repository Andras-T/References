package hu.bme.aut.steptracker.fragments

import android.os.Bundle
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.fragment.app.Fragment
import hu.bme.aut.steptracker.databinding.FragmentAchievementBinding
import hu.bme.aut.steptracker.model.Score
import hu.bme.aut.steptracker.sqlite.PersistentDataHelper

class AchievementFragment : Fragment() {

    private lateinit var binding: FragmentAchievementBinding
    private lateinit var dataHelper: PersistentDataHelper

    override fun onCreateView(inflater: LayoutInflater, container: ViewGroup?, savedInstanceState: Bundle?): View {
        super.onCreateView(inflater, container, savedInstanceState)
        binding = FragmentAchievementBinding.inflate(layoutInflater, container, false)
        dataHelper = activity?.let { PersistentDataHelper(it) }!!
        return binding.root
    }

    override fun onViewCreated(view: View, savedInstanceState: Bundle?) {
        super.onViewCreated(view, savedInstanceState)
        loadData()
    }

    private fun loadData() {
        dataHelper.open()
        var steps = 0
        val scores = dataHelper.getScore()
        var max = 0f
        for (score: Score in scores){
            steps += score.steps.toInt()
            if (score.steps > max)
                max = score.steps
        }

        if (max > 6000) {
            binding.achievement1.isChecked=true
            if (max > 12000) {
                binding.achievement2.isChecked=true
                if (max > 18000) {
                    binding.achievement3.isChecked=true
                    if(max > 24000) {
                        binding.achievement4.isChecked=true
                    }
                }
            }
        }
        dataHelper.close()
    }
}